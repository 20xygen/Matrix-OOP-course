#include <iostream>
#include "matrix.h"

template<int M, int N, typename Field>
void print(const Matrix<M, N, Field>& matrix) {
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      std::cout << matrix[i, j] << ' ';
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template<int M, int N, typename Field>
void print_int(const Matrix<M, N, Field>& matrix) {
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      std::cout << int(matrix[i, j]) << ' ';
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}



int main() {

  Residue<5> a(5);
  Residue<5> b(10);
  Residue<5> c = a / b;
  std::cout << int(c) << '\n';

  Matrix<4, 4, double> m0;
  auto m1 = Matrix<4, 4, double>::unityMatrix();
  m0 = m1 * 4.;
  print<4, 4, double>(m0);
  print<4, 4, double>(m0.inverted());
  std::cout << m0.det() << '\n';
  m0.invert();
  print<4, 4, double>(m0);

  Matrix<5, 5, double> m2 = Matrix<5, 5, double>::unityMatrix() * 5.;
  Matrix<5, 5, double> m3 = Matrix<5, 5, double>::unityMatrix() * 2.;
  print<5, 5, double>(m3);
  m3.invert();
  Matrix<5, 5, double> m4 = m2 * m3;
  print<5, 5, double>(m2);
  print<5, 5, double>(m3);
  print<5, 5, double>(m4);
  std::cout << m0.rank() << '\n';
  std::cout << m4.rank() << '\n';
  Matrix<4, 4, double> m5;
  std::cout << m5.rank() << '\n';

  const Matrix<2, 3, Residue<6>> m6 = {{1, 2, 3}, {0, 2, 4}};
  print_int<2, 3, Residue<6>>(m6);

  std::cout << "\n\n";

  Residue<99527> a3(17095);
  Residue<99527> a2(53);
  Residue<99527> a1 = a2 * a3;
  std::cout << int(a1) << '\n';
  std::cout << (17095 * 53) % 99527 << '\n';

  Residue<99527> r1(10292);
  Residue<99527> r2(53);
  Residue<99527> r3 = r1 / r2;
  std::cout << int(r3) << '\n';
  std::cout << (8568 * 53) % 99527 << '\n';
  std::cout << int(Residue<99527>(1) / r2) << '\n';
  std::cout << int((Residue<99527>(1) / r2) * r2) << '\n';
  std::cout << (62832 * 53) % 99527 << '\n';
  std::cout << int(r3 * r2) << '\n';

  Residue<13> q1 = 5;
  Residue<13> q2 = 8;
  Residue<13> q3 = q1 / q2;
  std::cout << int(q3) << '\n';
  std::cout << (8 * 12) % 13 << '\n';


  Matrix<2, 3, Residue<6>> w1 = {{10, 2, 3}, {0, 2, 4}};
  std::cout << (w1[1, 2] == Residue<6>(4)) << '\n';
  Matrix<3, 3, Residue<6>> w2 = Matrix<3, 3, Residue<6>>::unityMatrix();
  auto w3 = w2 * Residue<6>(4);
  std::cout << (w3[0, 0] == Residue<6>(4)) << '\n';
  w3[0, 0] = Residue<6>(3);
  std::cout << (w3[0, 0] == Residue<6>(4)) << '\n';
  std::cout << (w3[0, 0] == Residue<6>(3)) << '\n';

  Matrix<2, 2, Residue<7>> w4 = {{1, 2}, {0, 5}};
  Matrix<2, 3, Residue<7>> w5 = {{4, 2, 3}, {1, 2, 4}};
  auto w6 = w4 * w5;
  print_int<2, 3, Residue<7>>(w6);
  Matrix<4, 4, double> w7 = {{0, 4, 10, 1}, {4, 8, 18, 7}, {10, 18, 20, 17}, {1, 7, 17, 3}};
  std::cout << w7.rank() << '\n';
  Matrix<3, 5, double> w8 = {{2, -1, 3, -2, 4}, {4, -2, 5, 1, 7}, {2, -1, 1, 8, 2}};
  std::cout << w8.rank() << '\n';
  Matrix<5, 4, double> w9 = {{4, -7, -2, 1}, {-1, 3, 3, -4}, {-3, 5, 1, 0}, {-2, 3, 0, 1}, {1, -2, -1, 1}};
  std::cout << w9.rank() << '\n';
  Matrix<3, 4, double> w10 = {{1, 2, 3, 0}, {4, 5, 6, 0}, {7, 8, 9, 0}};
  std::cout << w10.rank() << '\n';
  std::cout << Matrix<10, 10, double>::unityMatrix().rank() << '\n';
  Matrix<3, 3, Residue<7>> w11 = {{1, 2, 3}, {1, 2, 0}, {4, 1, 5}};
  std::cout << w11.rank() << '\n';
  Matrix<3, 3, Residue<7>> w12;
  std::cout << w12.rank() << '\n';
  Matrix<5, 4, Residue<17>> w13 = {{4, 0, 3, 2}, {1, 11, 4, 5}, {7, 1, 5, 3}, {13, 15, 15, 0}, {1, 13, 2, 3}};
  std::cout << w13.rank() << '\n';

  Matrix<2, 2, double> w16 = {{1, 2}, {4, 5}};
  Matrix<2, 2, double> w17 = w16.inverted();
  print<2, 2, double>(w16);
  print<2, 2, double>(w17);

  Matrix<3, 3, Rational> w14 = {{1, 2, 3}, {1, 2, 0}, {4, 1, 5}};
  Matrix<3, 3, Rational> w15 = w14.inverted();
  print<3, 3, Rational>(w14);
  print<3, 3, Rational>(w15);

  Matrix<3, 3, double> w18 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 10}};
  Matrix<3, 3, double> w19 = w18.invertedLU();
  print<3, 3, double>(w19);

  Matrix<3, 3, Rational> w20 = {{1, 2, 3}, {1, 2, 0}, {4, 1, 5}};
  Matrix<3, 3, Rational> w21 = w20.invertedLU();
  print<3, 3, Rational>(w21);

  // test_inverse();

  return 0;
}
