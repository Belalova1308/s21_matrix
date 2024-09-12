#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}
S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 1 || cols < 1) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  matrix_ = new double[rows_ * cols_]();
}
S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  matrix_ = new double[other.rows_ * other.cols_]();
  std::copy(other.matrix_, other.matrix_ + rows_ * cols_, matrix_);
}
S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), cols_(other.cols_) {
  matrix_ = new double[other.rows_ * other.cols_]();
  std::copy(other.matrix_, other.matrix_ + rows_ * cols_, matrix_);
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}
S21Matrix::~S21Matrix() {
  if (matrix_) {
    delete[] matrix_;
  }
}
int S21Matrix::get_rows() const { return rows_; }
int S21Matrix::get_cols() const { return cols_; }
void S21Matrix::set_rows(const int &new_rows) {
  if (new_rows <= 0) throw std::length_error("Array size can't be zero");

  S21Matrix tmp(new_rows, cols_);
  for (int i = 0; i < (rows_ < new_rows ? rows_ : new_rows); ++i)
    for (int j = 0; j < cols_; ++j) tmp[i][j] = (*this)[i][j];

  *this = std::move(tmp);
}
void S21Matrix::set_cols(const int &new_cols) {
  if (new_cols <= 0) throw std::length_error("Array size can't be zero");

  S21Matrix tmp(rows_, new_cols);
  for (int i = 0; i < rows_; ++i)
    for (int j = 0; j < (cols_ < new_cols ? cols_ : new_cols); ++j)
      tmp[i][j] = (*this)[i][j];

  *this = std::move(tmp);
}
S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this != &other) {
    delete[] matrix_;

    rows_ = other.rows_;
    cols_ = other.cols_;

    matrix_ = new double[rows_ * cols_]();
    std::copy(other.matrix_, other.matrix_ + rows_ * cols_, matrix_);
  }
  return *this;
}
double *S21Matrix::operator[](int row) const {
  if (row >= rows_)
    throw std::out_of_range("Incorrect input, index is out of range");
  return &matrix_[row * cols_];
}
double &S21Matrix::operator()(int row, int col) {
  if (row >= rows_ || col >= cols_ || row < 0 || col < 0) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  return matrix_[row * cols_ + col];
}
double &S21Matrix::operator()(int row, int col) const {
  if (row >= rows_ || col >= cols_ || row < 0 || col < 0) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  return matrix_[row * cols_ + col];
}
S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix tmp{*this};
  tmp.SumMatrix(other);
  return tmp;
}
S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}
S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix tmp{*this};
  tmp.SubMatrix(other);
  return tmp;
}
S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}
S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix tmp{*this};
  tmp.MulMatrix(other);
  return tmp;
}
S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}
S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix tmp{*this};
  tmp.MulNumber(num);
  return tmp;
}
S21Matrix &S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}
bool S21Matrix::operator==(const S21Matrix &other) const {
  return EqMatrix(other);
}
void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        (*this)[i][j] += other[i][j];
      }
    }
  }
}
void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        (*this)[i][j] -= other[i][j];
      }
    }
  }
}
void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      (*this)[i][j] *= num;
    }
  }
}
void S21Matrix::MulMatrix(const S21Matrix &other) {
  S21Matrix res((*this).rows_, other.cols_);
  if ((*this).cols_ != other.rows_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have a.cols=b.rows");
  } else {
    for (int i = 0; i < res.rows_; i++) {
      for (int j = 0; j < res.cols_; j++) {
        double temp = 0.0;
        for (int k = 0; k < other.rows_; k++) {
          temp += (*this)[i][k] * other[k][j];
        }
        res[i][j] += temp;
      }
    }
  }
  *this = std::move(res);
}
S21Matrix S21Matrix::Transpose() {
  S21Matrix res((*this).cols_, (*this).rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res[j][i] = (*this)[i][j];
    }
  }
  return res;
}
double S21Matrix::Determinant() {
  double res = 0;
  if (rows_ != cols_) {
    throw std::out_of_range("Incorrect input, matrix is not square");
  } else if (rows_ == 1) {
    res = (*this)[0][0];
  } else if (rows_ == 2) {
    res = (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
  } else if (rows_ > 2) {
    int sing = 1;
    for (int i = 0; i < cols_; i++) {
      S21Matrix minor = CreateMinor(*this, 0, i);
      double minor_determinant = minor.Determinant();
      res += sing * (*this)[0][i] * minor_determinant;
      sing = -sing;
    }
  }
  return res;
}

S21Matrix S21Matrix::CreateMinor(const S21Matrix &m, int row_idx, int col_idx) {
  S21Matrix res{m.rows_ - 1, m.cols_ - 1};
  int minor_rows = 0;
  for (int i = 0; i < m.rows_; i++) {
    if (i == row_idx) continue;
    int minor_columns = 0;
    for (int j = 0; j < m.cols_; j++) {
      if (j == col_idx) continue;
      res[minor_rows][minor_columns] = m[i][j];
      minor_columns++;
    }
    minor_rows++;
  }
  return res;
}
S21Matrix S21Matrix::CalcComplements() {
  S21Matrix res{rows_, cols_};
  res.rows_ = this->rows_;
  res.cols_ = this->cols_;
  if (rows_ != cols_) {
    throw std::out_of_range("Incorrect input, matrices isn't square");
  } else {
    for (int i = 0; i < rows_; i++) {
      int sign = i % 2 == 0 ? 1 : -1;
      for (int j = 0; j < cols_; j++) {
        S21Matrix minor = CreateMinor(*this, i, j);
        double determinant = minor.Determinant();
        res[i][j] = sign * determinant;
        sign = -sign;
      }
    }
  }
  return res;
}
bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (std::abs(other(i, j) - (*this)(i, j)) > 1e-7) {
        return false;
      }
    }
  }
  return true;
}
S21Matrix S21Matrix::InverseMatrix() {
  S21Matrix res{rows_, cols_};
  res.rows_ = this->rows_;
  res.cols_ = this->cols_;
  double determinant = this->Determinant();
  if (determinant == 0) {
    throw std::out_of_range("Determinant is zero");
  } else {
    if (rows_ == 1) {
      res[0][0] = 1 / determinant;
    } else {
      S21Matrix matrixMinors = this->CalcComplements();
      S21Matrix matrixMinorsTransposed = matrixMinors.Transpose();
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < rows_; j++) {
          res[i][j] = matrixMinorsTransposed[i][j] * (1 / determinant);
        }
      }
    }
  }
  return res;
}