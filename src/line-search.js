'use strict';

/**
 * Fits a parabola with three points.
 * @ignore
 * @param {number} a1 The left value.
 * @param {number} a2 The middle value.
 * @param {number} a3 The right value.
 * @param {number} f1 The function value at a1.
 * @param {number} f2 The function value at a2.
 * @param {number} f3 The function value at a3.
 * @returns {number} The minimized x-axis point.
 */
function fitParabola (a1, a2, a3, f1, f2, f3) {
  const a23 = a2 - a3;
  const a31 = a3 - a1;
  const a12 = a1 - a2;
  const a23sq = a2 * a2 - a3 * a3;
  const a31sq = a3 * a3 - a1 * a1;
  const a12sq = a1 * a1 - a2 * a2;

  const den = 2 * (f1 * a23 + f2 * a31 + f3 * a12);
  if (den === 0) {
    return undefined;
  } else {
    const num = f1 * a23sq + f2 * a31sq + f3 * a12sq;
    return num / den;
  }
}

/**
 * Uses a parabolic fit to determine a minimizing step. Step delta is doubled each iteration until
 * the minimum is passed, then the last point is taken between the last two points, making the
 * steps equal. Then the parabola is fit and the minimum is calculated.
 * @ignore
 * @param {Vector2|Vector3|Vector4} X0 The starting point.
 * @param {Vector2|Vector3|Vector4} s The search direction.
 * @param {Function} F The objective function.
 * @param {number} delta The beginning numerical step.
 * @param {Vector2|Vector3|Vector4} XMin The placeholder for the minimized point.
 * @returns {Vector2|Vector3|Vector4} The minimized point.
 */
function parabolicLineSearch (X0, S, F, delta = 2.2e-16, XMin) {
  const s = S.clone().normalize();
  const X = XMin || X0.clone();
  const alphas = [0, 0, 0];
  const fs = [0, 0, 0];
  let j = 0;
  let alpha = delta;

  // find 3 points that bracket the minimum
  while (j < 3 || fs[(j - 2) % 3] - fs[(j - 1) % 3] > 0) {
    X.copy(X0).addScaledVector(s, alpha);
    fs[j % 3] = F(X);
    alphas[j % 3] = alpha;
    alpha *= 2;
    j++;
  }
  let a1 = alphas[j % 3];
  let f1 = fs[j % 3];
  let a2 = alphas[(j + 1) % 3];
  let f2 = fs[(j + 1) % 3];
  let a3 = alphas[(j - 1) % 3];
  let f3 = fs[(j - 1) % 3];

  // close in on the minimum to see if the parabolic fit still holds
  let aMin0 = Number.POSITIVE_INFINITY;
  let aMin1 = fitParabola(a1, a2, a3, f1, f2, f3);
  let MIN_TOL = 1e-8;
  let minFitIters = 0;
  let MIN_FIT_ITER_MAX = 2;
  // try a few times to shrink the boundaries and refit parabola
  while (minFitIters < MIN_FIT_ITER_MAX || Math.abs(aMin1 - aMin0) > MIN_TOL) {
    minFitIters++;
    // move the outermost boundary closer to the minimum
    if (f3 - f2 > f1 - f2) {
      // move the right boundary
      let ministep = 1;
      let aNew, fNew;
      do {
        // get progressively closer to the original a3 while not losing the minimum
        ministep++;
        aNew = a2 + ((a3 - a2) * (ministep - 1) / ministep);
        fNew = F(X.copy(X0).addScaledVector(s, aNew));
      } while (f2 > fNew);
      a3 = aNew;
      f3 = fNew;
    } else {
      // move the left boundary
      let ministep = 1;
      let aNew, fNew;
      do {
        // get progressively closer to the original a1 while not losing the minimum
        ministep++;
        aNew = a2 + ((a1 - a2) * (ministep - 1) / ministep);
        fNew = F(X.copy(X0).addScaledVector(s, aNew));
      } while (f2 > fNew);
      a1 = aNew;
      f1 = fNew;
    }
    aMin0 = aMin1;
    aMin1 = fitParabola(a1, a2, a3, f1, f2, f3);
  }
  return X.copy(X0).addScaledVector(s, aMin1);
}

module.exports = parabolicLineSearch;
