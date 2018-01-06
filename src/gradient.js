'use strict';

function hypot (a, b) {
  if (a === 0 && b === 0) {
    return 0;
  }
  var x = Math.abs(a);
  var y = Math.abs(b);
  var t = Math.min(x, y);
  var u = Math.max(x, y);
  t = t / u;
  return u * Math.sqrt(1 + t * t);
}

/**
 * Computes the central difference gradient. Returns the gradient vector 2-norm.
 * @ignore
 * @param {(Vector2|Vector3|Vector4)} X0 The evaluation point.
 * @param {number} h The delta in X0 to compute partial derivatives.
 * @param {Function} F The function to take the derivative of.
 * @param {(Vector2|Vector3|Vector4)} G Optional. A vector placeholder for the gradient.
 * @returns {number} The gradient vector 2-norm.
 */
function getGradient (X0, h, F, G) {
  let fx0 = 0;
  let fxh = 0;
  let xi = 0;
  const dim = X0.dimension;

  let nrm2 = 0;
  let gradVal = 0;
  for (let i = 0; i < dim; ++i) {
    xi = X0.getComponent(i);
    X0.setComponent(i, xi - h);
    fx0 = F(X0);
    X0.setComponent(i, xi + h);
    fxh = F(X0);
    gradVal = (fxh - fx0) / (2 * h);
    nrm2 = hypot(nrm2, gradVal);
    G.setComponent(i, gradVal);
  }
  if (isNaN(nrm2)) {
    throw new Error('2-norm of gradient is NaN.');
  }
  return nrm2;
}

module.exports = getGradient;
