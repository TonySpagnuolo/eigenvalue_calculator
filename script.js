const calculateButton = document.getElementById("calculate");
const results = document.getElementById("results");
const errorMessage = document.getElementById("error-message");
const normalizeButton = document.getElementById("normalize");

const grid = document.getElementById("matrix-grid");
const sizeLabel = document.getElementById("size-label");
let n = 0;

let eigenvectors = [];   // one entry per eigenvalue: an array of complex numbers, or null
const SUB = "₁₂₃₄";

document.getElementById("shrink").addEventListener("click", () => resize(n - 1));
document.getElementById("grow").addEventListener("click", () => resize(n + 1));

function resize(newN) {
    if (newN < 2 || newN > 4) return;   // clamp to 2-4

    // Remember what the user already typed so it survives the rebuild.
    const old = [];
    for (let i = 0; i < n; i++) {
        old.push([]);
        for (let j = 0; j < n; j++) {
            old[i].push(document.getElementById(`a-${i}-${j}`).value);
        }
    }

    n = newN;
    grid.style.gridTemplateColumns = `repeat(${n}, auto)`;
    grid.innerHTML = "";

    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            const input = document.createElement("input");
            input.type = "number";
            input.step = "any";
            input.id = `a-${i}-${j}`;
            input.value = old[i]?.[j] ?? "0";   // "?." handles the new rows/columns
            grid.appendChild(input);
        }
    }
    eigenvectors = [];

    sizeLabel.textContent = `${n}×${n}`;
    results.hidden = true;
    errorMessage.textContent = "";
}

function matMul(A, B) {
    const size = A.length;
    const C = A.map(() => new Array(size).fill(0));
    for (let i = 0; i < size; i++)
        for (let j = 0; j < size; j++)
            for (let k = 0; k < size; k++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

function trace(M) {
    let sum = 0;
    for (let i = 0; i < M.length; i++) sum += M[i][i];
    return sum;
}

// Returns [c0, c1, ..., cn] for p(λ) = c0 + c1·λ + ... + cn·λⁿ
function charPoly(A) {
    const size = A.length;
    const c = new Array(size + 1).fill(0);
    c[size] = 1;
    let M = A.map(row => row.map(() => 0));       // M₀ = 0

    for (let k = 1; k <= size; k++) {
        M = matMul(A, M);                          // Mₖ = A·Mₖ₋₁ + cₙ₋ₖ₊₁·I
        for (let i = 0; i < size; i++) M[i][i] += c[size - k + 1];
        c[size - k] = -trace(matMul(A, M)) / k;    // cₙ₋ₖ = −tr(A·Mₖ)/k
    }
    return c;
}

function readMatrix() {
    const A = [];
    for (let i = 0; i < n; i++) {
        A.push([]);
        for (let j = 0; j < n; j++) {
            const text = document.getElementById(`a-${i}-${j}`).value.trim();
            const value = Number(text);
            if (text === "" || !Number.isFinite(value)) {
                errorMessage.textContent = "Please enter a valid number in every matrix entry.";
                results.hidden = true;
                return null;
            }
            A[i].push(value);
        }
    }
    errorMessage.textContent = "";
    return A;
}

let normalized = false;

calculateButton.addEventListener("click", calculateEigenvalues);
normalizeButton.addEventListener("click", normalizeEigenvectors);



// Complex arithmetic on [re, im] pairs
const cAdd = (a, b) => [a[0] + b[0], a[1] + b[1]];
const cSub = (a, b) => [a[0] - b[0], a[1] - b[1]];
const cMul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const cDiv = (a, b) => {
    const d = b[0] * b[0] + b[1] * b[1];
    return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
};
const cAbs = (a) => Math.hypot(a[0], a[1]);

// Evaluate p(z) = c0 + c1·z + ... + cn·zⁿ using Horner's method
function polyEval(c, z) {
    let p = [c[c.length - 1], 0];
    for (let k = c.length - 2; k >= 0; k--) {
        p = cAdd(cMul(p, z), [c[k], 0]);
    }
    return p;
}

// Durand-Kerner: finds all n roots of a monic polynomial at once.
function findRoots(c) {
    const deg = c.length - 1;

    // Starting guesses: points spiraling out from the origin.
    // They must NOT be real, or they can never leave the real axis.
    let r = 1;
    for (let i = 0; i < deg; i++) r = Math.max(r, 1 + Math.abs(c[i]));
    const roots = [];
    let z = [r, 0];
    for (let k = 0; k < deg; k++) {
        roots.push(z);
        z = cMul(z, [0.4, 0.9]);
    }

    for (let iter = 0; iter < 500; iter++) {
        let maxStep = 0;
        for (let i = 0; i < deg; i++) {
            // denominator = product of (zᵢ − zⱼ) for all j ≠ i
            let denom = [1, 0];
            for (let j = 0; j < deg; j++) {
                if (j !== i) denom = cMul(denom, cSub(roots[i], roots[j]));
            }
            const step = cDiv(polyEval(c, roots[i]), denom);
            roots[i] = cSub(roots[i], step);
            maxStep = Math.max(maxStep, cAbs(step));
        }
        if (maxStep < 1e-12) break;
    }

    // Clean up: tiny imaginary parts are numerical noise.
    return roots.map(([re, im]) =>
        Math.abs(im) < 1e-7 * (1 + Math.abs(re)) ? [re, 0] : [re, im]
    );
}

// Merge roots within `tol` of each other into [value, multiplicity] pairs.
function groupRoots(roots, tol = 1e-5) {
    const groups = [];
    for (const r of roots) {
        const g = groups.find(g => cAbs(cSub(g.sum.map(x => x / g.count), r)) < tol);
        if (g) { g.sum = cAdd(g.sum, r); g.count++; }
        else groups.push({ sum: r, count: 1 });
    }
    return groups.map(g => ({ value: [g.sum[0] / g.count, g.sum[1] / g.count], mult: g.count }));
}

// Returns a list of basis vectors (each an array of [re, im]) for the null space of A − λI.
function nullSpace(A, lambda) {
    const size = A.length;
    const M = A.map((row, i) => row.map((x, j) => (i === j ? cSub([x, 0], lambda) : [x, 0])));

    let norm = 0;
    for (const row of A) for (const x of row) norm = Math.max(norm, Math.abs(x));
    const tol = 1e-6 * Math.max(norm, 1);

    // Gauss-Jordan elimination with partial pivoting.
    const pivotCols = [];
    let row = 0;
    for (let col = 0; col < size && row < size; col++) {
        let best = row;
        for (let r = row + 1; r < size; r++) {
            if (cAbs(M[r][col]) > cAbs(M[best][col])) best = r;
        }
        if (cAbs(M[best][col]) < tol) continue;          // no pivot: free column
        [M[row], M[best]] = [M[best], M[row]];

        const p = M[row][col];
        M[row] = M[row].map(x => cDiv(x, p));            // scale pivot to 1
        for (let r = 0; r < size; r++) {
            if (r === row) continue;
            const f = M[r][col];
            M[r] = M[r].map((x, j) => cSub(x, cMul(f, M[row][j])));
        }
        pivotCols.push(col);
        row++;
    }

    // One basis vector per free column.
    const basis = [];
    for (let free = 0; free < size; free++) {
        if (pivotCols.includes(free)) continue;
        const v = Array.from({ length: size }, () => [0, 0]);
        v[free] = [1, 0];
        pivotCols.forEach((pc, r) => { v[pc] = [-M[r][free][0], -M[r][free][1]]; });
        basis.push(v);
    }
    return basis;
}

// Prints a complex number like 3, 1.5 + 2i, -i, dropping zero parts.
function formatComplex(z) {
    const re = Number(z[0].toFixed(6));
    const im = Number(z[1].toFixed(6));
    if (im === 0) return `${re}`;
    const imText = Math.abs(im) === 1 ? "i" : `${Math.abs(im)}i`;
    if (re === 0) return im < 0 ? `-${imText}` : imText;
    return `${re} ${im < 0 ? "-" : "+"} ${imText}`;
}

// Scale so the first non-negligible entry is 1 (the "original" display form).
function scaleToFirstOne(v) {
    const p = v.find(z => cAbs(z) > 1e-6);
    return p ? v.map(z => cDiv(z, p)) : v;
}

// Scale to length 1: sqrt(sum of |entry|²).
function normalizeVector(v) {
    if (v === null) return null;
    const mag = Math.sqrt(v.reduce((s, z) => s + z[0] * z[0] + z[1] * z[1], 0));
    return mag === 0 ? null : v.map(z => [z[0] / mag, z[1] / mag]);
}

// Returns one value and one vector (or null) per eigenvalue, repeats included.
function solveEigen(A) {
    const groups = groupRoots(findRoots(charPoly(A)));
    groups.sort((x, y) => y.value[0] - x.value[0] || y.value[1] - x.value[1]);

    const values = [];
    const vectors = [];
    for (const g of groups) {
        let basis;
        if (g.value[1] < 0) {
            // Conjugate partner: reuse its vectors, conjugated.
            const conj = [g.value[0], -g.value[1]];
            const partner = groups.find(h => cAbs(cSub(h.value, conj)) < 1e-5);
            basis = partner
                ? nullSpace(A, partner.value).map(v => v.map(z => [z[0], -z[1]]))
                : nullSpace(A, g.value);
        } else {
            basis = nullSpace(A, g.value);
        }
        for (let k = 0; k < g.mult; k++) {
            values.push(g.value);
            vectors.push(basis[k] ? scaleToFirstOne(basis[k]) : null);
        }
    }
    return { values, vectors };
}

function formatVector(v) {
    if (v === null) return "No independent eigenvector";
    return `<span class="eigenvector">${v.map(z => `[ ${formatComplex(z)} ]`).join("<br>")}</span>`;
}

function showEigenvectors() {
    document.getElementById("eigenvectors").innerHTML = eigenvectors.map((v, i) => `
        <p>
            v${SUB[i]} =
            <br>
            ${formatVector(normalized ? normalizeVector(v) : v)}
        </p>
    `).join("");
}

function calculateEigenvalues() {
    const A = readMatrix();
    if (A === null) return;

    const { values, vectors } = solveEigen(A);
    eigenvectors = vectors;
    normalized = false;
    normalizeButton.textContent = "Normalize Eigenvectors";

    document.getElementById("eigenvalues").innerHTML =
        values.map((z, i) => `<p>λ${SUB[i]} = ${formatComplex(z)}</p>`).join("");

    showEigenvectors();
    results.hidden = false;
}

function normalizeEigenvectors() {
    if (eigenvectors.length === 0) {
        errorMessage.textContent = "Please calculate the eigenvalues and eigenvectors first.";
        return;
    }
    normalized = !normalized;
    normalizeButton.textContent = normalized ? "Show Original Eigenvectors" : "Normalize Eigenvectors";
    showEigenvectors();
}

resize(2);