/**
 * solver.js — Numerical Linear Algebra Engine
 * Eigenvalues via QR algorithm (power iteration for small matrices)
 * Eigenvectors via inverse iteration
 * Characteristic polynomial via Faddeev-LeVerrier algorithm
 */

"use strict";

const EigenSolver = (() => {

  /* ── Helpers ── */

  const eps = 1e-10;

  function matMul(A, B) {
    const n = A.length;
    const C = zeros(n);
    for (let i = 0; i < n; i++)
      for (let k = 0; k < n; k++)
        for (let j = 0; j < n; j++)
          C[i][j] += A[i][k] * B[k][j];
    return C;
  }

  function matScale(A, s) {
    return A.map(row => row.map(v => v * s));
  }

  function matAdd(A, B) {
    return A.map((row, i) => row.map((v, j) => v + B[i][j]));
  }

  function matSub(A, B) {
    return A.map((row, i) => row.map((v, j) => v - B[i][j]));
  }

  function eye(n) {
    return Array.from({length: n}, (_, i) =>
      Array.from({length: n}, (_, j) => i === j ? 1 : 0));
  }

  function zeros(n) {
    return Array.from({length: n}, () => new Array(n).fill(0));
  }

  function trace(A) {
    return A.reduce((s, row, i) => s + row[i], 0);
  }

  function matCopy(A) {
    return A.map(row => [...row]);
  }

  function vecNorm(v) {
    return Math.sqrt(v.reduce((s, x) => s + x * x, 0));
  }

  function vecNormalize(v) {
    const n = vecNorm(v);
    return n < eps ? v : v.map(x => x / n);
  }

  function vecScale(v, s) { return v.map(x => x * s); }
  function vecDot(a, b) { return a.reduce((s, x, i) => s + x * b[i], 0); }
  function vecSub(a, b) { return a.map((x, i) => x - b[i]); }

  /* ── Determinant (LU) ── */
  function det(A) {
    const n = A.length;
    const L = matCopy(A);
    let sign = 1;
    for (let i = 0; i < n; i++) {
      let maxRow = i;
      for (let k = i + 1; k < n; k++)
        if (Math.abs(L[k][i]) > Math.abs(L[maxRow][i])) maxRow = k;
      if (maxRow !== i) { [L[i], L[maxRow]] = [L[maxRow], L[i]]; sign *= -1; }
      if (Math.abs(L[i][i]) < eps) return 0;
      for (let k = i + 1; k < n; k++) {
        const f = L[k][i] / L[i][i];
        for (let j = i; j < n; j++) L[k][j] -= f * L[i][j];
      }
    }
    return sign * L.reduce((p, row, i) => p * row[i], 1);
  }

  /* ── Rank ── */
  function rank(A) {
    const n = A.length;
    const M = matCopy(A);
    let r = 0;
    for (let col = 0; col < n && r < n; col++) {
      let pivot = -1;
      for (let row = r; row < n; row++)
        if (Math.abs(M[row][col]) > eps) { pivot = row; break; }
      if (pivot === -1) continue;
      [M[r], M[pivot]] = [M[pivot], M[r]];
      const scale = M[r][col];
      for (let j = col; j < n; j++) M[r][j] /= scale;
      for (let row = 0; row < n; row++) {
        if (row === r) continue;
        const f = M[row][col];
        for (let j = col; j < n; j++) M[row][j] -= f * M[r][j];
      }
      r++;
    }
    return r;
  }

  /* ── Characteristic Polynomial via Faddeev-LeVerrier ── */
  function characteristicPoly(A) {
    const n = A.length;
    const coeffs = new Array(n + 1).fill(0);
    coeffs[n] = 1;
    let M = eye(n);
    for (let k = 1; k <= n; k++) {
      M = matAdd(matMul(A, M), matScale(eye(n), coeffs[n - k + 1]));
      coeffs[n - k] = -trace(M) / k;
    }
    return coeffs; // coeffs[0]*λ^0 + ... + coeffs[n]*λ^n
  }

  /* ── Complex Number Helpers ── */
  function cAdd(a, b) { return { r: a.r + b.r, i: a.i + b.i }; }
  function cSub(a, b) { return { r: a.r - b.r, i: a.i - b.i }; }
  function cMul(a, b) { return { r: a.r*b.r - a.i*b.i, i: a.r*b.i + a.i*b.r }; }
  function cDiv(a, b) {
    const d = b.r*b.r + b.i*b.i;
    return { r: (a.r*b.r + a.i*b.i)/d, i: (a.i*b.r - a.r*b.i)/d };
  }
  function cAbs(a) { return Math.sqrt(a.r*a.r + a.i*a.i); }
  function cSqrt(a) {
    const m = Math.sqrt(cAbs(a));
    const theta = Math.atan2(a.i, a.r) / 2;
    return { r: m * Math.cos(theta), i: m * Math.sin(theta) };
  }
  function cCbrt(a) {
    const m = Math.cbrt(cAbs(a));
    const theta = Math.atan2(a.i, a.r) / 3;
    return { r: m * Math.cos(theta), i: m * Math.sin(theta) };
  }
  function c(r, i = 0) { return { r, i }; }

  /* ── Eigenvalues ── */

  // 2×2 analytic
  function eig2(A) {
    const a = A[0][0], b = A[0][1], cc = A[1][0], d = A[1][1];
    const tr = a + d;
    const disc = tr*tr - 4*(a*d - b*cc);
    if (disc >= 0) {
      const sq = Math.sqrt(disc);
      return [{ r: (tr + sq)/2, i: 0 }, { r: (tr - sq)/2, i: 0 }];
    }
    const sq = Math.sqrt(-disc);
    return [{ r: tr/2, i: sq/2 }, { r: tr/2, i: -sq/2 }];
  }

  // Cubic roots (Cardano)
  function solveCubic(a, b, cc, d) {
    // Normalize: x^3 + px^2 + qx + r = 0
    const p = b/a, q = cc/a, r = d/a;
    const p2 = p/3, p3 = p2*p2*p2;
    const Q = (3*q - p*p)/9;
    const R = (9*p*q - 27*r - 2*p*p*p)/54;
    const D = Q*Q*Q + R*R;
    const roots = [];
    if (D > 1e-10) {
      const sqrtD = Math.sqrt(D);
      const S = Math.cbrt(R + sqrtD);
      const T = Math.cbrt(R - sqrtD);
      roots.push({ r: S + T - p2, i: 0 });
      roots.push({ r: -(S+T)/2 - p2, i: Math.sqrt(3)*(S-T)/2 });
      roots.push({ r: -(S+T)/2 - p2, i: -Math.sqrt(3)*(S-T)/2 });
    } else if (D < -1e-10) {
      const theta = Math.acos(R / Math.sqrt(-Q*Q*Q));
      const sqQ = 2 * Math.sqrt(-Q);
      roots.push({ r: sqQ*Math.cos(theta/3) - p2, i: 0 });
      roots.push({ r: sqQ*Math.cos((theta + 2*Math.PI)/3) - p2, i: 0 });
      roots.push({ r: sqQ*Math.cos((theta + 4*Math.PI)/3) - p2, i: 0 });
    } else {
      const S = Math.cbrt(R);
      roots.push({ r: 2*S - p2, i: 0 });
      roots.push({ r: -S - p2, i: 0 });
      roots.push({ r: -S - p2, i: 0 });
    }
    return roots;
  }

  // 3×3 analytic
  function eig3(A) {
    const p = characteristicPoly(A);
    // p[0] + p[1]λ + p[2]λ² + p[3]λ³ = 0  → reverse
    return solveCubic(p[3], p[2], p[1], p[0]);
  }

  // QR iteration for n >= 4
  function eigQR(A) {
    const n = A.length;
    let H = matCopy(A);
    const maxIter = 2000;

    for (let iter = 0; iter < maxIter; iter++) {
      // Wilkinson shift
      const a = H[n-2][n-2], b = H[n-2][n-1];
      const cc2 = H[n-1][n-2], d = H[n-1][n-1];
      const tr = a + d;
      const delta = (a - d) / 2;
      const sign = delta >= 0 ? 1 : -1;
      const mu = d - sign * (b*cc2) / (Math.abs(delta) + Math.sqrt(delta*delta + b*cc2));

      // Shifted QR step
      let Q = eye(n);
      let R = matCopy(H);
      for (let k = 0; k < n; k++) {
        for (let l = n - 1; l > k; l--) {
          if (Math.abs(R[l][k]) < eps) continue;
          const r = Math.hypot(R[l-1][k], R[l][k]);
          const cos = R[l-1][k] / r;
          const sin = -R[l][k] / r;
          // Apply Givens rotation
          for (let j = 0; j < n; j++) {
            const t1 = R[l-1][j], t2 = R[l][j];
            R[l-1][j] = cos*t1 - sin*t2;
            R[l][j] = sin*t1 + cos*t2;
          }
          for (let j = 0; j < n; j++) {
            const t1 = Q[j][l-1], t2 = Q[j][l];
            Q[j][l-1] = cos*t1 - sin*t2;
            Q[j][l] = sin*t1 + cos*t2;
          }
        }
      }
      H = matMul(R, Q);
    }

    // Extract eigenvalues from quasi-upper-triangular H
    const eigenvalues = [];
    let i = 0;
    while (i < n) {
      if (i < n - 1 && Math.abs(H[i+1][i]) > eps) {
        const [e1, e2] = eig2([[H[i][i], H[i][i+1]], [H[i+1][i], H[i+1][i+1]]]);
        eigenvalues.push(e1, e2);
        i += 2;
      } else {
        eigenvalues.push({ r: H[i][i], i: 0 });
        i++;
      }
    }
    return eigenvalues;
  }

  /* ── Eigenvectors via Inverse Iteration ── */
  function eigenvector(A, lambda) {
    const n = A.length;
    // Shift: (A - λI)
    const shift = lambda.r;
    const Ashifted = matCopy(A);
    for (let i = 0; i < n; i++) Ashifted[i][i] -= (shift + 1e-7);

    // Random initial vector
    let v = Array.from({length: n}, (_, i) => i === 0 ? 1 : 0.1 * (i+1));
    v = vecNormalize(v);

    // Power iteration on (A-λI)^{-1} → solve (A-λI)x = v repeatedly
    for (let iter = 0; iter < 80; iter++) {
      const next = luSolve(Ashifted, v);
      if (!next) break;
      const newNorm = vecNorm(next);
      if (newNorm < eps) break;
      v = next.map(x => x / newNorm);
    }
    return v;
  }

  /* ── LU Solve ── */
  function luSolve(A, b) {
    const n = A.length;
    const M = A.map((row, i) => [...row, b[i]]);
    for (let i = 0; i < n; i++) {
      let maxRow = i;
      for (let k = i + 1; k < n; k++)
        if (Math.abs(M[k][i]) > Math.abs(M[maxRow][i])) maxRow = k;
      [M[i], M[maxRow]] = [M[maxRow], M[i]];
      if (Math.abs(M[i][i]) < eps) return null;
      for (let k = i + 1; k < n; k++) {
        const f = M[k][i] / M[i][i];
        for (let j = i; j <= n; j++) M[k][j] -= f * M[i][j];
      }
    }
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      x[i] = M[i][n];
      for (let j = i + 1; j < n; j++) x[i] -= M[i][j] * x[j];
      x[i] /= M[i][i];
    }
    return x;
  }

  /* ── Format Helpers ── */
  function formatNum(x, digits = 5) {
    if (Math.abs(x) < 1e-9) return '0';
    if (Number.isInteger(x) || Math.abs(x - Math.round(x)) < 1e-8)
      return String(Math.round(x));
    return parseFloat(x.toFixed(digits)).toString();
  }

  function formatComplex(c, digits = 5) {
    const r = Math.abs(c.r) < 1e-9 ? 0 : c.r;
    const im = Math.abs(c.i) < 1e-9 ? 0 : c.i;
    if (im === 0) return formatNum(r, digits);
    const imStr = (Math.abs(im) === 1) ? '' : formatNum(Math.abs(im), digits);
    const sign = im < 0 ? ' − ' : ' + ';
    if (r === 0) return (im < 0 ? '−' : '') + imStr + 'i';
    return formatNum(r, digits) + sign + imStr + 'i';
  }

  function formatPoly(coeffs) {
    const n = coeffs.length - 1;
    let parts = [];
    for (let k = n; k >= 0; k--) {
      const c = coeffs[k];
      const rounded = Math.round(c * 1e8) / 1e8;
      if (Math.abs(rounded) < 1e-7) continue;
      const absC = Math.abs(rounded);
      const numStr = Math.abs(absC - 1) < 1e-9 && k > 0 ? '' : formatNum(absC);
      let term = '';
      if (k === 0) term = numStr || '1';
      else if (k === 1) term = (numStr ? numStr : '') + 'λ';
      else term = (numStr ? numStr : '') + 'λ' + toSup(k);
      parts.push({ sign: rounded < 0 ? '−' : '+', term, k });
    }
    if (parts.length === 0) return '0';
    let result = '';
    parts.forEach((p, idx) => {
      if (idx === 0) {
        result += (p.sign === '−' ? '−' : '') + p.term;
      } else {
        result += ' ' + p.sign + ' ' + p.term;
      }
    });
    return result + ' = 0';
  }

  function toSup(n) {
    const sups = '⁰¹²³⁴⁵';
    return String(n).split('').map(d => sups[d] || d).join('');
  }

  /* ── Main Solve Entry Point ── */
  function solve(matrix) {
    const n = matrix.length;
    let eigenvalues;
    if (n === 2) eigenvalues = eig2(matrix);
    else if (n === 3) eigenvalues = eig3(matrix);
    else eigenvalues = eigQR(matrix);

    // Deduplicate complex conjugates for display
    const evecs = eigenvalues.map(ev => {
      if (Math.abs(ev.i) > 1e-6) {
        // Complex eigenvector — return null, note it
        return null;
      }
      return eigenvector(matrix, ev);
    });

    const poly = characteristicPoly(matrix);
    const d = det(matrix);
    const tr = trace(matrix);
    const rk = rank(matrix);

    return {
      eigenvalues: eigenvalues.map(ev => formatComplex(ev)),
      eigenvectors: evecs.map((v, i) => ({
        label: 'v' + toSub(i+1),
        eigenvalue: formatComplex(eigenvalues[i]),
        isComplex: Math.abs(eigenvalues[i].i) > 1e-6,
        vector: v ? v.map(x => formatNum(x)) : null
      })),
      poly: formatPoly(poly),
      trace: formatNum(tr),
      det: formatNum(d),
      rank: rk
    };
  }

  function toSub(n) {
    const subs = '₀₁₂₃₄₅';
    return String(n).split('').map(d => subs[d] || d).join('');
  }

  return { solve, formatNum };

})();