'use strict';

/**
 * Updates the matrix with a BFGS rank-2 update.
 * @ignore
 * @param {(Matrix2|Matrix3|Matrix4)} N The matrix to update.
 * @param {(Vector2|Vector3|Vector4)} y The gradient difference vector.
 * @param {(Vector2|Vector3|Vector4)} dx The x vector difference.
 */
function bfgsUpdate (N, y, dx) {
  if (N.dimension !== y.dimension ||
      N.dimension !== dx.dimension ||
      y.dimension !== dx.dimension) {
    throw new Error('bfgsUpdate(): Dimension mismatch.');
  }
  var t1 = y.clone().multiplyMatrix(N);
  var a = dx.dot(y);
  var b = y.dot(t1);
  var c = 1.0 / a;
  var d = (1 + b / a) * c;
  N.addOuterProduct(dx, dx, d);
  N.addOuterProduct(dx, t1, -c);
  N.addOuterProduct(t1, dx, -c);
  return N;
}

module.exports = bfgsUpdate;
