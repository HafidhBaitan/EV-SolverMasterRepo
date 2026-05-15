# λ EigenSolve

> An interactive, browser-based eigenvalue and eigenvector calculator for matrices from 2×2 up to 5×5. No dependencies, no build tools — just open and solve.



---

## Features

- **Multi-size support** — switch between 2×2, 3×3, 4×4, and 5×5 matrices with a single click
- **Exact & numerical methods** — analytic solutions for 2×2 and 3×3; QR algorithm with Wilkinson shift for 4×4 and 5×5
- **Full output** — eigenvalues, eigenvectors, characteristic polynomial, trace, determinant, and rank
- **Complex eigenvalue detection** — complex eigenvalues are displayed in `a + bi` form; complex eigenvectors are flagged
- **Dark / Light mode** — toggle with preference saved to `localStorage`
- **Matrix presets** — Random, Identity, and Clear buttons for quick setup
- **Copy to clipboard** — export all results as plain text in one click
- **Keyboard friendly** — `Tab` navigates cells, `Enter` or `Ctrl+Enter` triggers solve
- **Fully responsive** — works on desktop, tablet, and mobile down to 320px wide
- **Zero dependencies** — pure HTML, CSS, and vanilla JavaScript; no frameworks or libraries

---

## Project Structure

```
eigensolve/
├── index.html   # Page structure and layout
├── styles.css   # All styling — design tokens, dark/light themes, animations
├── solver.js    # Numerical math engine (eigenvalue algorithms)
└── app.js       # UI controller — rendering, events, interactivity
```

---

## Getting Started

No installation or build step needed.

1. Download or clone all four files into the same folder.
2. Open `index.html` in any modern browser.
3. Enter matrix values, pick a size, and click **Solve**.

```bash
# Or serve locally with Python
python -m http.server 8080
# Then visit http://localhost:8080
```

---

## How It Works

### Algorithms by Matrix Size

| Size | Eigenvalue Method | Notes |
|------|-------------------|-------|
| 2×2  | Analytic (quadratic formula) | Exact; handles complex roots |
| 3×3  | Analytic (Cardano's cubic formula) | Exact for three roots including complex conjugate pairs |
| 4×4  | QR algorithm + Wilkinson shift | Iterative; converges to quasi-upper-triangular form |
| 5×5  | QR algorithm + Wilkinson shift | Same as 4×4; up to 2000 iterations |

### Characteristic Polynomial

Computed via the **Faddeev–LeVerrier algorithm**, which builds the polynomial coefficients iteratively using the matrix trace at each step — stable and works for any n×n.

### Eigenvectors

Computed via **inverse iteration** (power iteration on the shifted inverse matrix). For a given eigenvalue λ, the method repeatedly solves `(A − λI)x = v` and normalizes, converging to the dominant eigenvector.

### Additional Properties

- **Determinant** — LU decomposition with partial pivoting
- **Rank** — Gaussian elimination with tolerance `ε = 1e-10`
- **Trace** — direct diagonal sum (equals sum of eigenvalues)

---

## Keyboard Shortcuts

| Key | Action |
|-----|--------|
| `Tab` / `Shift+Tab` | Move between matrix cells |
| `Enter` (in a cell) | Trigger solve |
| `Ctrl+Enter` / `Cmd+Enter` | Trigger solve from anywhere |

---

## Browser Support

Works in all modern browsers (Chrome, Firefox, Safari, Edge). Requires ES6+ support. No polyfills included.

---

## Limitations

- Matrix entries must be real numbers (complex input is not supported)
- For 4×4 and 5×5 matrices, eigenvalues are numerical approximations — results are accurate to ~5 significant figures
- Degenerate matrices (e.g. all-zero rows) may produce `0` or `NaN` for some values
- Complex eigenvectors are detected but not numerically computed — only the eigenvalue is shown for those cases

---

## Design

Built with a scientific / terminal aesthetic using:

- **Syne** (display font) — geometric, high-contrast headings
- **JetBrains Mono** — monospaced entries, values, and results
- CSS custom properties for seamless dark/light theming
- Animated bracket notation around the matrix input grid
- Staggered cell entrance animations on size change

---

## License

MIT — free to use, modify, and distribute.

---

*EigenSolve — powered by numerical linear algebra, no libraries required.*
