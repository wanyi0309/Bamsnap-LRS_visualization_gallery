#!/usr/bin/env python3
"""
Dump a compact per-read intermediate state from BamSnap-LRS immediately before
DNA read rendering.

This script does NOT re-parse or re-implement BamSnap-LRS alignment handling.
It directly calls BamSnap-LRS's existing reader/layout/renderer helpers to
capture the information needed to inspect:
  - alignment identity/status
  - query-space grouping and ordering
  - stack-row assignment
  - read-level color/opacity/orientation
  - compact segment-level rendering semantics

The TSV contains one row per visible BAM alignment record (Read object).
Only compact rendering semantics are stored as JSON; pixel geometry and full
Read/Segment object dumps are intentionally omitted.
"""

import argparse
import csv
import json
import os
import sys

# ---------------------------------------------------------------------------
# Bootstrap import path BEFORE importing bamsnap_lrs.
#
# Bamsnap-LRS uses a src/ layout:
#   /path/to/Bamsnap-LRS/
#       src/
#           bamsnap_lrs/
#
# Therefore --bamsnap-src must point to the directory that CONTAINS
# the bamsnap_lrs package, i.e. ".../Bamsnap-LRS/src".
# ---------------------------------------------------------------------------
_bootstrap = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
_bootstrap.add_argument("--bamsnap-src", default=None)
_bootstrap_args, _ = _bootstrap.parse_known_args()

if _bootstrap_args.bamsnap_src:
    src_dir = os.path.abspath(_bootstrap_args.bamsnap_src)
    if not os.path.isdir(src_dir):
        raise SystemExit(f"--bamsnap-src does not exist or is not a directory: {src_dir}")
    sys.path.insert(0, src_dir)

try:
    from bamsnap_lrs.reader import fetch_reads
except ModuleNotFoundError as e:
    if e.name == "bamsnap_lrs":
        raise SystemExit(
            "Cannot import bamsnap_lrs.\n"
            "Please provide --bamsnap-src pointing to Bamsnap-LRS/src, e.g.\n"
            "  --bamsnap-src /data/work/01.bamsnap_lrs/Bamsnap-LRS/src"
        ) from e
    raise

from bamsnap_lrs.layout import (
    assign_split_read_stacks,
    query_order_key,
    _query_span_from_segments,
    segments_to_pixels,
)
from bamsnap_lrs.svg_renderer import (
    _build_supplementary_read_styles,
    _merge_pixel_spans,
    rgb_to_hex,
    MATCH_HEX,
    AUTO_SIMPLIFY_BP_PER_PX,
    AUTO_MIN_DEL_BP,
    AUTO_MIN_DEL_PX,
)
from bamsnap_lrs.styles import (
    color_for_type,
    STRAND_COLORS,
    shade_by_mapq,
    MISMATCH_COLORS,
)


# Match svg_renderer.render_svg_snapshot() default left/right margin.
MARGIN = 20


def parse_pos(text):
    chrom, coords = text.split(":", 1)
    coords = coords.replace(",", "")
    if "-" in coords:
        s, e = coords.split("-", 1)
        return chrom, int(s), int(e)
    pos = int(coords)
    return chrom, max(0, pos - 250), pos + 250


def jd(obj):
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"))


def base_color_for_rect(rect_type, read, color_by):
    if color_by == "type":
        rgb = color_for_type(rect_type)
    elif color_by == "strand":
        rgb = STRAND_COLORS["-" if read.reverse else "+"]
    elif color_by == "mapq":
        rgb = shade_by_mapq(color_for_type(rect_type), read.mapq)
    else:
        # "base" is handled explicitly for mismatch segments below; other
        # segment types fall back to their type color, matching the renderer.
        rgb = color_for_type(rect_type)
    return rgb_to_hex(rgb)


def build_render_semantics(
    read,
    rects,
    margin,
    width,
    bp_per_px,
    detail,
    color_by,
    supp_style,
    auto_simplify,
    hide_insertions,
    show_insertion_labels,
):
    """
    Reproduce the non-highlight DNA read-body decisions made in svg_renderer,
    but return only compact semantic information rather than pixel geometry.

    Each returned item contains:
      type, length, symbol, color, opacity, label, drawn
    """
    supp_color = rgb_to_hex(supp_style[0]) if supp_style else None
    supp_opacity = float(supp_style[1]) if supp_style else None
    default_opacity = supp_opacity if supp_opacity is not None else 1.0

    # Map renderer rect index back to parsed Segment index exactly as the
    # renderer does. Insertions need special handling because they consume no
    # reference bases.
    rect_to_seg = {}
    rect_idx = 0
    for seg_idx, seg in enumerate(read.segments):
        if seg.ref_consumed == 0:
            if seg.type == "ins":
                if rect_idx < len(rects) and rects[rect_idx][0] == "ins":
                    rect_to_seg[rect_idx] = seg_idx
                    rect_idx += 1
            continue
        if rect_idx < len(rects):
            rect_to_seg[rect_idx] = seg_idx
            rect_idx += 1

    # At overview scale, supplementary match/mismatch/small-deletion pieces
    # are merged into one or more body spans. We keep only the semantic fact
    # that those merged bodies are drawn; their x-coordinates are deliberately
    # omitted from this compact output.
    merged_rect_indices = set()
    merged_body_count = 0
    if supp_style and auto_simplify:
        body_spans = []
        for ridx, (render_type, x0, x1) in enumerate(rects):
            x0d = max(margin, min(width - margin, margin + x0))
            x1d = max(margin, min(width - margin, margin + x1))
            if x1d <= x0d:
                continue

            body_like = render_type in ("match", "mismatch")
            if render_type == "del":
                seg_idx = rect_to_seg.get(ridx)
                seg = read.segments[seg_idx] if seg_idx is not None else None
                if seg is not None:
                    del_bp = int(getattr(seg, "length", seg.ref_consumed))
                    body_like = (
                        del_bp < AUTO_MIN_DEL_BP
                        or (x1d - x0d) < AUTO_MIN_DEL_PX
                    )

            if body_like:
                body_spans.append((x0d, x1d))
                merged_rect_indices.add(ridx)

        merged_body_count = len(_merge_pixel_spans(body_spans))

    semantics = []

    # Record the actual merged supplementary body primitive(s) that remain
    # visible after overview simplification.
    for _ in range(merged_body_count):
        semantics.append({
            "type": "supplementary_overview_body",
            "length": None,
            "symbol": "rect",
            "color": supp_color,
            "opacity": default_opacity,
            "label": None,
            "drawn": True,
        })

    for ridx, (render_type, x0, x1) in enumerate(rects):
        seg_idx = rect_to_seg.get(ridx)
        seg = read.segments[seg_idx] if seg_idx is not None else None
        seg_length = getattr(seg, "length", None)

        x0d = max(margin, min(width - margin, margin + x0))
        x1d = max(margin, min(width - margin, margin + x1))

        item = {
            "type": render_type,
            "length": seg_length,
            "symbol": "rect",
            "color": None,
            "opacity": default_opacity,
            "label": None,
            "drawn": True,
        }

        # The individual primitive is not drawn separately if it has been
        # merged into the supplementary overview body recorded above.
        if ridx in merged_rect_indices:
            item["drawn"] = False
            item["color"] = supp_color
            semantics.append(item)
            continue

        color_hex = base_color_for_rect(render_type, read, color_by)

        if render_type == "ins":
            if auto_simplify:
                item["drawn"] = False
                item["symbol"] = None
                semantics.append(item)
                continue

            item["color"] = MATCH_HEX if hide_insertions else color_hex
            if (
                not hide_insertions
                and detail == "high"
                and show_insertion_labels
                and seg is not None
                and seg.length > 0
            ):
                item["label"] = f"I({seg.length})"

        elif render_type == "ref_skip":
            item["symbol"] = "line"
            item["color"] = "#b0c4de"

        elif render_type == "del":
            small_overview_del = False
            if auto_simplify and seg is not None:
                del_bp = int(getattr(seg, "length", seg.ref_consumed))
                del_width_px = max(0, x1d - x0d)
                small_overview_del = (
                    del_bp < AUTO_MIN_DEL_BP
                    or del_width_px < AUTO_MIN_DEL_PX
                )

            if small_overview_del:
                item["color"] = supp_color if supp_style else MATCH_HEX
            else:
                item["color"] = "#808080"
                if (
                    detail == "high"
                    and show_insertion_labels
                    and seg is not None
                    and seg.length > 0
                ):
                    item["label"] = str(seg.length)

        elif render_type == "mismatch":
            if auto_simplify:
                item["color"] = supp_color if supp_style else MATCH_HEX
            else:
                # A merged mismatch segment is colored using its first query
                # base, matching the existing renderer behavior.
                mismatch_color = MATCH_HEX
                if read.seq and seg_idx is not None:
                    read_cursor = sum(
                        s.read_consumed for s in read.segments[:seg_idx]
                    )
                    if 0 <= read_cursor < len(read.seq):
                        base = read.seq[read_cursor].upper()
                        mismatch_color = rgb_to_hex(
                            MISMATCH_COLORS.get(base, (200, 60, 60))
                        )
                item["color"] = mismatch_color

        elif render_type == "match":
            item["color"] = supp_color if supp_style else MATCH_HEX

        else:
            if supp_style and render_type in ("match", "soft", "hard"):
                item["color"] = supp_color
            else:
                item["color"] = color_hex

        semantics.append(item)

    return semantics


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Dump a compact BamSnap-LRS per-read intermediate state immediately "
            "before DNA read rendering."
        )
    )
    ap.add_argument(
        "--bamsnap-src",
        default=None,
        help="Path to Bamsnap-LRS/src. If omitted, bamsnap_lrs must already be importable.",
    )
    ap.add_argument("--bam", required=True, help="Input BAM/CRAM used for plotting")
    ap.add_argument("--pos", required=True, help="Region, e.g. chr14:54808983-54821699")
    ap.add_argument("--out", required=True, help="Output TSV path")

    # Match the relevant BamSnap-LRS DNA defaults/filtering behavior.
    ap.add_argument("--mapq", type=int, default=0)
    ap.add_argument("--show-supp", action="store_true")
    ap.add_argument("--show-secondary", action="store_true")
    ap.add_argument("--use-md", action="store_true")
    ap.add_argument("--use-cs", action="store_true")
    ap.add_argument("--fa", default=None)

    # These rendering arguments are still needed internally because
    # render_semantics_json should reflect the same rendering decisions as
    # BamSnap-LRS. They are intentionally not repeated as TSV columns.
    ap.add_argument("--width", type=int, default=1200)
    ap.add_argument("--detail", choices=["low", "mid", "high"], default="mid")
    ap.add_argument(
        "--overview-detail",
        choices=["hide", "show"],
        default="hide",
    )
    ap.add_argument(
        "--color-by",
        choices=["type", "base", "strand", "mapq"],
        default="type",
    )
    ap.add_argument("--hide-insertions", action="store_true")
    ap.add_argument(
        "--no-insertion-labels",
        dest="show_insertion_labels",
        action="store_false",
    )
    ap.set_defaults(show_insertion_labels=True)

    args = ap.parse_args()
    chrom, start, end = parse_pos(args.pos)

    reads = fetch_reads(
        args.bam,
        chrom,
        start,
        end,
        mapq_min=args.mapq,
        show_supp=args.show_supp,
        show_secondary=args.show_secondary,
        use_md=args.use_md,
        use_cs=args.use_cs,
        use_ref=bool(args.fa),
        fa_path=args.fa,
    )

    if not reads:
        raise SystemExit(f"No reads returned for {chrom}:{start}-{end}")

    # Exact split-read stack assignment used by the DNA rendering path.
    stacks = assign_split_read_stacks(reads, start, end)

    # Group visible alignment records by qname, matching the current
    # BamSnap-LRS supplementary styling/query-order path.
    groups = {}
    for idx, read in enumerate(reads):
        groups.setdefault(read.qname, []).append(idx)

    supplementary_styles = _build_supplementary_read_styles(
        reads, groups, stacks
    )

    # Query-order rank within each visible qname group.
    query_rank = {}
    group_size = {}
    for qname, idxs in groups.items():
        ordered = sorted(idxs, key=lambda i: query_order_key(reads[i]))
        for rank, idx in enumerate(ordered, start=1):
            query_rank[idx] = rank
            group_size[idx] = len(ordered)

    content_width = args.width - 2 * MARGIN
    if content_width <= 0:
        raise SystemExit(f"--width must be greater than {2 * MARGIN}")
    bp_per_px = (end - start) / float(content_width)

    auto_simplify = (
        bp_per_px > AUTO_SIMPLIFY_BP_PER_PX
        and args.overview_detail != "show"
    )

    fields = [
        "qname",
        "chrom",
        "start",
        "end",
        "strand",
        "mapq",
        "primary",
        "supplementary",
        "secondary",
        "query_start",
        "query_end",
        "query_length",
        "layout_query_start",
        "layout_query_end",
        "qname_group_size",
        "query_order_in_group",
        "stack_row",
        "supplementary_group_color",
        "supplementary_group_opacity",
        "read_body_color",
        "read_body_opacity",
        "arrow_direction",
        "render_semantics_json",
    ]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for idx, read in enumerate(reads):
            supp_style = supplementary_styles.get(idx)
            supp_color = rgb_to_hex(supp_style[0]) if supp_style else None
            supp_opacity = float(supp_style[1]) if supp_style else None

            layout_qs, layout_qe = _query_span_from_segments(read)

            rects = segments_to_pixels(
                read.segments,
                read.start,
                start,
                bp_per_px,
                detail=args.detail,
            )

            render_semantics = build_render_semantics(
                read=read,
                rects=rects,
                margin=MARGIN,
                width=args.width,
                bp_per_px=bp_per_px,
                detail=args.detail,
                color_by=args.color_by,
                supp_style=supp_style,
                auto_simplify=auto_simplify,
                hide_insertions=args.hide_insertions,
                show_insertion_labels=args.show_insertion_labels,
            )

            # The visible match/read body uses supplementary group styling when
            # present; otherwise the current SVG renderer uses MATCH_HEX.
            body_color = supp_color if supp_style else MATCH_HEX
            body_opacity = supp_opacity if supp_style else 1.0

            writer.writerow({
                "qname": read.qname,
                "chrom": read.chrom,
                "start": read.start,
                "end": read.end,
                "strand": "-" if read.reverse else "+",
                "mapq": read.mapq,
                "primary": read.primary,
                "supplementary": read.supplementary,
                "secondary": read.secondary,
                "query_start": read.query_start,
                "query_end": read.query_end,
                "query_length": read.query_length,
                "layout_query_start": layout_qs,
                "layout_query_end": layout_qe,
                "qname_group_size": group_size[idx],
                "query_order_in_group": query_rank[idx],
                "stack_row": stacks[idx],
                "supplementary_group_color": supp_color or ".",
                "supplementary_group_opacity": (
                    f"{supp_opacity:.6f}" if supp_opacity is not None else "."
                ),
                "read_body_color": body_color,
                "read_body_opacity": f"{body_opacity:.6f}",
                "arrow_direction": "left" if read.reverse else "right",
                "render_semantics_json": jd(render_semantics),
            })

    print(f"Wrote {len(reads)} alignment records to {args.out}")


if __name__ == "__main__":
    main()
