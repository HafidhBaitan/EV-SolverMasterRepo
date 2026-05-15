/**
 * app.js — EigenSolve UI Controller
 */

"use strict";

(function () {

  /* ── State ── */
  let currentSize = 3;

  /* ── DOM References ── */
  const matrixGrid   = document.getElementById('matrixGrid');
  const solvBtn      = document.getElementById('solveBtn');
  const clearBtn     = document.getElementById('clearBtn');
  const randomBtn    = document.getElementById('randomBtn');
  const identityBtn  = document.getElementById('identityBtn');
  const copyBtn      = document.getElementById('copyBtn');
  const themeToggle  = document.getElementById('themeToggle');
  const resultsSection = document.getElementById('resultsSection');
  const polyDisplay  = document.getElementById('polyDisplay');
  const eigenvaluesGrid = document.getElementById('eigenvaluesGrid');
  const eigenvectorsGrid = document.getElementById('eigenvectorsGrid');
  const traceVal     = document.getElementById('traceVal');
  const detVal       = document.getElementById('detVal');
  const rankVal      = document.getElementById('rankVal');
  const toast        = document.getElementById('toast');
  const sizeTabs     = document.querySelectorAll('.size-tab');

  /* ── Theme ── */
  const savedTheme = localStorage.getItem('eigensolve-theme') || 'dark';
  setTheme(savedTheme);

  themeToggle.addEventListener('click', () => {
    const curr = document.documentElement.getAttribute('data-theme');
    setTheme(curr === 'dark' ? 'light' : 'dark');
  });

  function setTheme(theme) {
    document.documentElement.setAttribute('data-theme', theme);
    localStorage.setItem('eigensolve-theme', theme);
  }

  /* ── Size Tabs ── */
  sizeTabs.forEach(tab => {
    tab.addEventListener('click', () => {
      const size = parseInt(tab.dataset.size, 10);
      sizeTabs.forEach(t => t.classList.remove('active'));
      tab.classList.add('active');
      currentSize = size;
      buildMatrix(size);
      hideResults();
    });
  });

  /* ── Build Matrix Grid ── */
  function buildMatrix(n, values = null) {
    matrixGrid.innerHTML = '';
    matrixGrid.style.gridTemplateColumns = `repeat(${n}, var(--cell-size))`;

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const input = document.createElement('input');
        input.type = 'number';
        input.inputMode = 'decimal';
        input.step = 'any';
        input.className = 'cell-input highlight';
        input.setAttribute('aria-label', `Row ${i+1}, Column ${j+1}`);
        input.dataset.row = i;
        input.dataset.col = j;

        if (values && values[i] && values[i][j] !== undefined) {
          input.value = values[i][j];
        }

        // Stagger animation delay
        input.style.animationDelay = `${(i * n + j) * 18}ms`;

        // Tab navigation override: natural row-by-row
        input.addEventListener('keydown', (e) => {
          if (e.key === 'Enter') { e.preventDefault(); solveMatrix(); }
          if (e.key === 'Tab') {
            const idx = i * n + j;
            const next = matrixGrid.querySelectorAll('.cell-input')[idx + (e.shiftKey ? -1 : 1)];
            if (next) { e.preventDefault(); next.focus(); }
          }
        });

        matrixGrid.appendChild(input);
      }
    }

    // Update bracket heights
    const cellSize = parseInt(getComputedStyle(document.documentElement).getPropertyValue('--cell-size')) || 56;
    const gap = n >= 5 ? 4 : 6;
    const bracketH = n * cellSize + (n - 1) * gap;
    document.documentElement.style.setProperty('--bracket-h', bracketH + 'px');
  }

  /* ── Get Matrix Values ── */
  function getMatrix() {
    const inputs = matrixGrid.querySelectorAll('.cell-input');
    const A = [];
    for (let i = 0; i < currentSize; i++) {
      A.push([]);
      for (let j = 0; j < currentSize; j++) {
        const val = parseFloat(inputs[i * currentSize + j].value);
        if (isNaN(val)) return null;
        A[i].push(val);
      }
    }
    return A;
  }

  /* ── Presets ── */
  clearBtn.addEventListener('click', () => {
    matrixGrid.querySelectorAll('.cell-input').forEach(inp => inp.value = '');
    hideResults();
  });

  randomBtn.addEventListener('click', () => {
    const n = currentSize;
    const vals = Array.from({length: n}, () =>
      Array.from({length: n}, () => Math.round((Math.random() * 20 - 10) * 10) / 10)
    );
    buildMatrix(n, vals);
    hideResults();
  });

  identityBtn.addEventListener('click', () => {
    const n = currentSize;
    const vals = Array.from({length: n}, (_, i) =>
      Array.from({length: n}, (_, j) => i === j ? 1 : 0)
    );
    buildMatrix(n, vals);
    hideResults();
  });

  /* ── Solve ── */
  solvBtn.addEventListener('click', solveMatrix);

  function solveMatrix() {
    const A = getMatrix();
    if (!A) {
      showToast('Please fill in all matrix entries.', true);
      shakeGrid();
      return;
    }

    // Button loading state
    solvBtn.classList.add('loading');
    solvBtn.querySelector('span:first-child').textContent = 'Solving…';

    // Defer for paint
    setTimeout(() => {
      try {
        const result = EigenSolver.solve(A);
        renderResults(result);
      } catch (err) {
        showToast('Computation error: ' + err.message, true);
        console.error(err);
      }
      solvBtn.classList.remove('loading');
      solvBtn.querySelector('span:first-child').textContent = 'Solve';
    }, 30);
  }

  /* ── Render Results ── */
  function renderResults(r) {
    // Polynomial
    polyDisplay.innerHTML = formatPolyHTML(r.poly);

    // Eigenvalues
    eigenvaluesGrid.innerHTML = '';
    r.eigenvalues.forEach((ev, i) => {
      const chip = document.createElement('div');
      chip.className = 'eigenvalue-chip';
      chip.style.animationDelay = `${i * 60}ms`;
      chip.innerHTML = `
        <span class="eigenvalue-label">λ${['₁','₂','₃','₄','₅'][i] || (i+1)}</span>
        <span class="eigenvalue-val">${escapeHTML(ev)}</span>
      `;
      eigenvaluesGrid.appendChild(chip);
    });

    // Eigenvectors
    eigenvectorsGrid.innerHTML = '';
    r.eigenvectors.forEach((ev, i) => {
      const block = document.createElement('div');
      block.className = 'eigenvector-block';
      block.style.animationDelay = `${i * 80}ms`;

      let innerHTML = `<div class="evec-label">For λ${['₁','₂','₃','₄','₅'][i] || (i+1)} = ${escapeHTML(ev.eigenvalue)}</div>`;

      if (ev.isComplex) {
        innerHTML += `<div class="evec-bracket-wrap">
          <div style="color:var(--text-dimmer);font-size:12px;padding:8px 0;">Complex eigenvector<br><span style="color:var(--text-dimmer);font-size:10px;">(use complex arithmetic)</span></div>
        </div>`;
      } else if (ev.vector) {
        const entries = ev.vector.map(x => `<div class="evec-entry">${escapeHTML(x)}</div>`).join('');
        innerHTML += `<div class="evec-bracket-wrap">
          <span class="evec-bracket">⎡</span>
          <div class="evec-column">${entries}</div>
          <span class="evec-bracket">⎤</span>
        </div>`;
      }

      block.innerHTML = innerHTML;
      eigenvectorsGrid.appendChild(block);
    });

    // Metrics
    traceVal.textContent = r.trace;
    detVal.textContent = r.det;
    rankVal.textContent = r.rank;

    // Show
    resultsSection.classList.add('visible');
    resultsSection.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  }

  function formatPolyHTML(poly) {
    // Highlight λ terms
    return poly
      .replace(/λ([⁰¹²³⁴⁵]?)/g, '<span class="poly-term-accent">λ$1</span>')
      .replace(/([\d.]+)(?=<span)/g, '<span class="poly-coeff">$1</span>');
  }

  function escapeHTML(s) {
    return String(s)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;');
  }

  /* ── Copy All ── */
  copyBtn.addEventListener('click', () => {
    const lines = [];
    lines.push('=== EigenSolve Results ===');
    lines.push('\nCharacteristic Polynomial:');
    lines.push(polyDisplay.textContent);
    lines.push('\nEigenvalues:');
    document.querySelectorAll('.eigenvalue-chip').forEach(chip => {
      lines.push('  ' + chip.textContent.trim());
    });
    lines.push('\nMatrix Properties:');
    lines.push(`  Trace: ${traceVal.textContent}`);
    lines.push(`  Determinant: ${detVal.textContent}`);
    lines.push(`  Rank: ${rankVal.textContent}`);

    navigator.clipboard.writeText(lines.join('\n'))
      .then(() => showToast('Results copied to clipboard ✓'))
      .catch(() => showToast('Copy failed — try manual selection', true));
  });

  /* ── Toast ── */
  let toastTimer;
  function showToast(msg, isError = false) {
    toast.textContent = msg;
    toast.className = 'toast show' + (isError ? ' error' : '');
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => toast.classList.remove('show'), 3000);
  }

  /* ── Shake Animation ── */
  function shakeGrid() {
    matrixGrid.animate([
      { transform: 'translateX(0)' },
      { transform: 'translateX(-6px)' },
      { transform: 'translateX(6px)' },
      { transform: 'translateX(-4px)' },
      { transform: 'translateX(4px)' },
      { transform: 'translateX(0)' }
    ], { duration: 300, easing: 'ease' });
  }

  /* ── Hide Results ── */
  function hideResults() {
    resultsSection.classList.remove('visible');
  }

  /* ── Keyboard shortcut: Ctrl+Enter to solve ── */
  document.addEventListener('keydown', (e) => {
    if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
      e.preventDefault();
      solveMatrix();
    }
  });

  /* ── Init ── */
  buildMatrix(currentSize);

  // Pre-fill a demo 3×3
  const demo = [[4, -2, 1], [2, 0, 1], [2, -2, 3]];
  setTimeout(() => buildMatrix(currentSize, demo), 100);

})();