'use strict';

const minimatrix = require('minimatrix');
const Matrix2 = minimatrix.Matrix2;
const Matrix3 = minimatrix.Matrix3;

const lineSearch = require('./line-search.js');
const gradient = require('./gradient.js');
const bfgsUpdate = require('./bfgsUpdate.js');

/**
 * The default values for solutions.
 * @ignore
 */
const defaults = Object.freeze({
  maxIterations: 200,
  tolerance: 1e-12
});

/**
 * @typedef {Object} ObjectiveOptions
 * @property {Vector2|Vector3|Vector4} start The start point.
 * @property {Function} func The objective function.
 * @property {number} delta The numerical step for the derivative.
 */

/**
 * @typedef {Object} SolutionOptions
 * @property {number} tolerance The numerical tolerance under which the function value is considered a solution.
 * @property {number} maxIterations The maximum number of iterations before the search is forcibly ended.
 */

/**
 * @typedef {Object} QuasiNewtonOptions
 * @property {ObjectiveOptions} objective The objective function options.
 * @property {SolutionOptions} solution The solution options.
 */

/**
 * @typedef OptimizationResults
 * @property {boolean} solutionValid True if the solution meets the criteria of an optimum, false if not.
 * @property {number} iterations The number of iterations taken.
 * @property {number} objective The value of the objective function.
 * @property {number} gradNorm The norm of the gradient at the solution.
 */

/**
 * @typedef ResultObject
 * @property {boolean} solutionValid True if the solution meets the criteria for a solution, false if not.
 * @property {number} iterations The number of iterations taken to arrive at the solution.
 * @property {number} objective The value of the objective function at the solution.
 * @property {number} gradNorm The norm of the gradient vector.
 */

/**
 * Performs an unconstrained quasi-Newton iteration on the objective function to find the minimum.
 * @memberof Optimization
 * @param {QuasiNewtonOptions} options Options for the optimization.
 * @returns {ResultObject} The results of the optimization.
 */
function quasiNewton (options) {
  if (!options.objective) {
    throw new Error('Undefined optimization objective.');
  }
  if (!options.objective.start) {
    throw new Error('Undefined start position.');
  }
  if (!options.objective.func) {
    throw new Error('Undefined objective function.');
  }

  var maxIterations = defaults.maxIterations;
  var tolerance = defaults.tolerance;
  if (options.solution) {
    if (options.solution.maxIterations &&
      !isNaN(options.solution.maxIterations)) {
      maxIterations = options.solution.maxIterations;
    } else {
      console.warn('Maximum iterations capped at default of ' + maxIterations + '.');
    }
    if (options.solution.tolerance &&
      !isNaN(options.solution.tolerance)) {
      tolerance = options.solution.tolerance;
    } else {
      console.warn('Numerical tolerance is default of ' + tolerance + '.');
    }
  }

  const F = options.objective.func;
  const G = (x, grad) => gradient(x, delta, F, grad);

  let x0 = options.objective.start.clone();
  const delta = options.objective.delta;
  const dim = x0.dimension;
  if (dim > 4 || dim < 2) {
    throw new Error('Dimension is out of range.');
  }
  let x1 = x0.clone();
  const dx = x1.clone();
  let grad0 = x0.clone();
  let grad1 = x0.clone();
  const y = x0.clone();
  const N = dim === 2 ? new Matrix2() : (dim === 3 ? new Matrix3() : undefined);
  let f0 = Number.POSITIVE_INFINITY;
  let f1 = F(x0);
  let gradNorm = G(x0, grad0);

  // use this for the search direction
  y.copy(grad0).multiplyScalar(-1.0 / gradNorm);

  let iter = 0;
  let temp;
  while (Math.abs(f1 - f0) > tolerance &&
         Math.abs(gradNorm) > tolerance &&
         iter++ <= maxIterations) {
    lineSearch(x0, y, F, 1e-15, x1);
    gradNorm = G(x1, grad1);
    f0 = f1;
    f1 = F(x1);
    dx.copy(x1).sub(x0);
    y.copy(grad1).sub(grad0);
    bfgsUpdate(N, y, dx);

    // get the new search direction
    y.copy(grad1).multiplyScalar(-1.0 / gradNorm).multiplyMatrix(N);
    // swap gradients
    temp = grad0;
    grad0 = grad1;
    grad1 = temp;
    // swap x values
    temp = x0;
    x0 = x1;
    x1 = temp;
  }
  return {
    solutionValid: Math.abs(f1 - f0) < tolerance || Math.abs(gradNorm) < tolerance,
    iterations: iter,
    solution: x0,
    objective: f1,
    gradNorm
  };
}

module.exports = quasiNewton;
