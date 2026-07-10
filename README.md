# Collaborating on the slides

This repository holds a slide deck built with [Quarto](https://quarto.org) using the
[reveal.js](https://quarto.org/docs/presentations/revealjs/) format. The slides are written in a
plain-text `.qmd` file (Quarto Markdown), and rendering that file produces a self-contained HTML
presentation you open in a browser.

> Throughout this README the source file is referred to as `slides.qmd`. Replace that with the
> actual file name in this repo.

## Before you start

1. Install Quarto: <https://quarto.org/docs/get-started/>
2. Get the repo (`git clone …`, or `git pull` if you already have it). **Always pull the latest
   changes before you start editing** so you're working on top of everyone else's work.

## Repository layout

```
slides.qmd          the presentation source (edit this)
css/nwpage.scss     the theme — imported in the YAML header, do not move
media/              images, logos, video, audio
data/               datasets and any other supporting files
```

A few conventions to keep the project tidy:

- **Theme:** the look of the deck comes from `css/nwpage.scss`, referenced in the slide header
  (see below). If you want to change fonts, colors, or spacing, edit that file rather than styling
  individual slides.
- **Media:** put every image, logo, or video in `media/` and reference it from there.
- **Everything else:** datasets, spreadsheets, and other supporting files go in `data/`.

## Add yourself to the authors

The deck's metadata lives in the YAML header at the very top of `slides.qmd`, between the two `---`
lines. Add your name to the `author` list:

```yaml
---
title: "Presentation title"
author:
  - name: Existing Author
  - name: Your Name          # <- add yourself here
format:
  revealjs:
    theme: css/nwpage.scss   # <- the theme; leave this pointing at css/nwpage.scss
---
```

Quarto uses [standard author metadata](https://quarto.org/docs/authoring/create-citeable-articles.html),
so each `- name:` entry can also carry optional fields like `email`, `affiliation`, or `orcid` if
you want them shown on the title slide.

## Creating a slide

Slides are separated by Markdown headings — you don't insert page breaks manually. See the
[presentation basics](https://quarto.org/docs/presentations/) docs for the full rundown.

- `##` (level-2 heading) starts a **new slide**, using the heading text as the slide title.
- `#` (level-1 heading) starts a **section / title slide** that introduces a group of slides.
- `---` on its own line inserts a slide break without a title.

```markdown
# Results                      <!-- section divider -->

## First finding               <!-- a slide -->

- Point one
- Point two

## Second finding              <!-- the next slide -->

Some text on this slide.
```

## Markdown basics

The slide body is [Markdown](https://quarto.org/docs/authoring/markdown-basics.html). The common
pieces:

```markdown
**bold text** and *italic text*

- bullet list
- another item
  - nested item

1. numbered list
2. second item

[a hyperlink](https://quarto.org)

> a block quote

`inline code`
```

For a block of code, wrap it in a line of three backticks above and below, with an optional
language name on the opening line (for example `r` or `python`). If you mark the block as a live
chunk with `{r}` or `{python}`, Quarto runs the code and drops the output onto the slide.

Images are just Markdown links to files in `media/`:

```markdown
![A short caption](media/example-plot.png)
```

To load data in a code chunk, point at the `data/` folder, e.g. `read_csv("data/results.csv")`.

## Using raw HTML

Anything Markdown can't express, you can write as raw HTML directly in the `.qmd` file and it passes
straight through to the slide:

```html
<div style="text-align: center; color: #0c5460;">
  <strong>Centered, colored text</strong>
</div>
```

Quarto/reveal.js also supports "fenced div" blocks for slide features (columns, speaker notes,
incremental lists, etc.). For example, side-by-side columns:

```markdown
::: {.columns}
::: {.column width="50%"}
Left side
:::
::: {.column width="50%"}
Right side
:::
:::
```

See the [reveal.js guide](https://quarto.org/docs/presentations/revealjs/) for the available
classes.

## Render from the command line

From the repo folder, build the HTML with:

```bash
quarto render slides.qmd
```

This produces `slides.html` alongside the source file. While actively editing, use live preview
instead — it opens the deck in your browser and reloads automatically every time you save:

```bash
quarto preview slides.qmd
```

If the project contains a `_quarto.yml` file, running `quarto render` with no file name builds the
whole project at once.

## View the slides

Open the rendered `slides.html` in any web browser. Double-click it in your file manager.

Use the arrow keys to navigate; press `S` for speaker notes, `F` for fullscreen, and `?` to see all
reveal.js keyboard shortcuts.

## Documentation

- Quarto home: <https://quarto.org>
- Presentations overview: <https://quarto.org/docs/presentations/>
- reveal.js format guide: <https://quarto.org/docs/presentations/revealjs/>
- Markdown basics: <https://quarto.org/docs/authoring/markdown-basics.html>
- reveal.js themes: <https://quarto.org/docs/presentations/revealjs/themes.html>
- Full reveal.js option reference: <https://quarto.org/docs/reference/formats/presentations/revealjs.html>