#ifndef __S21MATRIX_H__
#define __S21MATRIX_H__

#include <iostream>
class S21Matrix {
 private:
  int rows_, cols_;
  double *matrix_;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other);
  ~S21Matrix();
  double *operator[](int row) const;
  int get_cols() const;
  int get_rows() const;
  void set_rows(const int &new_rows);
  void set_cols(const int &new_cols);

  double &operator()(int row, int col);
  double &operator()(int row, int col) const;
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix operator+(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  bool operator==(const S21Matrix &other) const;
  S21Matrix operator*(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix operator*(const double num);
  S21Matrix &operator*=(const double num);
  S21Matrix &operator=(const S21Matrix &other);
  bool EqMatrix(const S21Matrix &other) const;
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix CreateMinor(const S21Matrix &m, int row_idx, int col_idx);
  S21Matrix InverseMatrix();
};

#endif
