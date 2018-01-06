/* global describe it */
'use strict';

const chai = require('chai');
const expect = chai.expect;
const minimatrix = require('minimatrix');
const Vector2 = minimatrix.Vector2;

const parabolicLineSearch = require('../src/line-search.js');
const quasiNewton = require('../src/quasiNewton.js');

describe('Basic Operations', () => {
  describe('Parabolic Line Search', () => {
    it('Test 1', () => {
      const TOL = 1e-15;
      const fParabola = function (X) {
        var x1 = X.x;
        var x2 = X.y;
        var f = x1 * x1 + x2 * x2;
        return f;
      };
      const X0 = new Vector2(-2, -2);
      const s = new Vector2(1, 1);
      const X = parabolicLineSearch(X0, s, fParabola, 0.01);
      expect(Math.abs(X.x)).to.be.below(TOL, 'X coordinate is not within tolerance.');
      expect(Math.abs(X.y)).to.be.below(TOL, 'Y coordinate is not within tolerance.');
    });
  });
});
describe('Optimization Problems', () => {
  it('Test 1', () => {
    const TOLERANCE = 1e-11;
    const MAX_ITERATIONS = 5;
    const F = function (X) {
      var x1 = X.x;
      var x2 = X.y;
      var f = x1 * x1 - 2 * x1 * x2 + 4 * x2 * x2;
      return f;
    };

    const x0 = new Vector2(-3, 1);
    const options = {
      objective: {
        start: x0,
        func: F,
        delta: 1e-13
      },
      solution: {
        tolerance: TOLERANCE,
        maxIterations: MAX_ITERATIONS
      }
    };
    const results = quasiNewton(options);
    chai.assert(results.solutionValid, 'Solution is not optimal');
    chai.assert.isBelow(results.objective, TOLERANCE, 'Result is not within tolerance');
    // chai.assert(results.iterations === 2, 'Should take only 2 iterations to reach optimum');
  });
});
