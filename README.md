# Bamsnap-LRS visual comparison gallery

A dependency-free static gallery for comparing matched long-read alignment visualizations produced by **Bamsnap-LRS**, **Wally**, and **SVhawkeye**. The gallery contains insertion, deletion, inversion, duplication, and RNA-splicing examples.

## Preview locally

You can open `index.html` directly in a modern browser. No build step or web framework is required.

## Publish with GitHub Pages

1. Add the contents of this directory to a GitHub repository.
2. Push the files to the default branch (for example, `main`).
3. In the repository, open **Settings → Pages**.
4. Under **Build and deployment**, choose **Deploy from a branch**.
5. Select the default branch and the **/(root)** folder, then save.
6. GitHub will provide the public Pages URL after deployment finishes.

Because all paths are relative, the site works both at a repository Pages subpath (`https://USER.github.io/REPO/`) and at a custom domain.

## Structure

- `index.html` — page markup
- `styles.css` — responsive page styling
- `app.js` — gallery interactions and deep linking
- `data.js` — matched locus/image manifest
- `assets/` — original visualization images, organized by class and tool
- `data/` — BED files supplied for the SV examples plus an RNA region index generated from image filenames

## Deep links

A specific example can be linked using query parameters, for example:

`?type=INS&id=INS1`

Supported `type` values are `INS`, `DEL`, `INV`, `DUP`, and `RNA`.

## Matching rules used to assemble the gallery

- INS and DEL: Bamsnap-LRS/SVhawkeye images are matched to the exact BED interval; Wally is matched by BED record ID (`INS1`, `DEL1`, etc.).
- INV and DUP: Bamsnap-LRS is matched to the displayed BED interval; SVhawkeye filenames encode the interval after removing the 2,000-bp flank from each side; Wally is matched by BED record ID.
- RNA: Bamsnap-LRS and Wally filenames encode the same displayed interval; SVhawkeye filenames encode that interval after removing the 100-bp flank from each side.

Missing outputs are intentionally retained as `Not available` rather than being replaced with a different locus. In the uploaded result set, Bamsnap-LRS images were absent for **INS27**, **INS42**, and **DUP28**.
