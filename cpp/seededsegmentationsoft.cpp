/****************************************************************************
**
** Copyright (C) 2020-2022 Harlen Batagelo <hbatagelo@gmail.com>,
**                         João Paulo Gois <jpgois@gmail.com>.
**
** This file is part of the implementation of the paper
** 'Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation'
** by Wallace Casaca, João Paulo Gois, Harlen Batagelo, Gabriel Taubin and
** Luis Gustavo Nonato. DOI 10.1109/TPAMI.2020.2974475.
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program. If not, see <http://www.gnu.org/licenses/>.
**
****************************************************************************/

#include "seededsegmentationsoft.h"

#include <Eigen/CholmodSupport>
#include <opencv2/opencv.hpp>

namespace lc {

SeededSegmentationSoft::SeededSegmentationSoft(QImage const *input,
                                               QImage const *seeds,
                                               QColor foreground,
                                               QColor background, double beta)
    : m_input{input}, m_cols{input->width()}, m_rows{input->height()},
      m_numPixels{m_rows * m_cols}, m_beta{beta} {
  Q_ASSERT(input != nullptr);
  Q_ASSERT(seeds != nullptr);

  // Split the different layers of the image and store them in three Eigen
  // matrices
  m_channels[0] = Eigen::MatrixXd(m_rows, m_cols);
  m_channels[1] = Eigen::MatrixXd(m_rows, m_cols);
  m_channels[2] = Eigen::MatrixXd(m_rows, m_cols);

  for (auto y{0}; y < m_rows; ++y) {
    for (auto x{0}; x < m_cols; ++x) {
      auto const color{input->pixelColor(x, y)};
      m_channels[0](y, x) = color.redF();
      m_channels[1](y, x) = color.greenF();
      m_channels[2](y, x) = color.blueF();
    }
  }

  // Set up seeds matrix
  m_seeds = Eigen::MatrixXd(m_rows, m_cols);

  for (auto y{0}; y < m_rows; ++y) {
    for (auto x{0}; x < m_cols; ++x) {
      auto const color{seeds->pixelColor(x, y)};
      if (color == foreground)
        m_seeds(y, x) = 1;
      else if (color == background)
        m_seeds(y, x) = -1;
      else
        m_seeds(y, x) = 0;
    }
  }
}

// Initialize the sparse matrices with the number of expected values in each row
// and set up the variables which will help us to compute the graph weight and
// the sum of the rows
void SeededSegmentationSoft::SetUpGraphWeightAndSum() {
  constexpr auto epsilon{1e-6};
  constexpr auto epsilonInv{1.0 / epsilon};

  std::array const neighborOffsets{
      Eigen::Vector2i{-1, 0}, Eigen::Vector2i{-1, -1}, Eigen::Vector2i{0, -1},
      Eigen::Vector2i{1, -1}, Eigen::Vector2i{1, 0},   Eigen::Vector2i{1, 1},
      Eigen::Vector2i{0, 1},  Eigen::Vector2i{-1, 1}};

  auto const n{m_numPixels};
  auto const valency{neighborOffsets.size() + 1};

  Eigen::Vector2i const boundsLo{0, 0};
  Eigen::Vector2i const boundsHi{m_rows, m_cols};

  m_sparseWeights = Eigen::SparseMatrix<double>(n, n);
  m_sparseWeights.reserve(Eigen::VectorXi::Constant(n, valency));

  m_sparseD = Eigen::SparseMatrix<double>(n, n);
  m_sparseD.reserve(Eigen::VectorXi::Constant(n, 1));

  // Compute the weight of each pixel with respect to its neighbors
  for (auto idx{0}; idx < n; ++idx) {
    auto summ{0.0};
    // Pixel row (x) and column (y)
    Eigen::Vector2i pI{idx / m_cols, idx % m_cols};

    for (auto const &offset : neighborOffsets) {
      Eigen::Vector2i const pJ{pI + offset};

      if (!inBounds(pJ, boundsLo, boundsHi))
        continue;

      // Compute difference of the values of each layer between the two pixels
      Eigen::Vector3d vecDiff{
          m_channels[0](pI.x(), pI.y()) - m_channels[0](pJ.x(), pJ.y()),
          m_channels[1](pI.x(), pI.y()) - m_channels[1](pJ.x(), pJ.y()),
          m_channels[2](pI.x(), pI.y()) - m_channels[2](pJ.x(), pJ.y())};

      // Calculate infinity norm of the differences
      auto const dist{vecDiff.lpNorm<Eigen::Infinity>()};
      // Compute the weight
      auto wij{std::exp(-m_beta * dist * dist)};
      wij = std::round(wij * epsilonInv) * epsilon;
      wij += epsilon;
      m_sparseWeights.insert(idx, (pJ.x() * m_cols) + pJ.y()) = wij;

      summ += wij;
    }

    // Insert the sum in the diagonal matrix
    m_sparseD.insert(idx, idx) = summ;
  }

  m_sparseWeights.makeCompressed();
  m_sparseD.makeCompressed();
}

// Set up the diagonal sparse matrix with ones if the pixel belongs to the
// background or the foreground, and zeros otherwise
void SeededSegmentationSoft::CheckIfInSeeds() {
  auto const n{m_numPixels};
  m_sparseIs = Eigen::SparseMatrix<double>(n, n);
  m_sparseIs.reserve(Eigen::VectorXi::Constant(n, 1));

  for (auto i{0}; i < m_rows; ++i) {
    for (auto j{0}; j < m_cols; ++j) {
      if (m_seeds(i, j) != 0.0) {
        auto const index{(i * m_cols) + j};
        m_sparseIs.insert(index, index) = 1;
      }
    }
  }

  m_sparseIs.makeCompressed();
}

// Compute difference between the row sum of the weights and the weights
void SeededSegmentationSoft::DiffBetweenSumAndWeights() {
  m_sparseL = m_sparseD;
  m_sparseL -= m_sparseWeights;
}

void SeededSegmentationSoft::Compute() {
  SetUpGraphWeightAndSum();
  CheckIfInSeeds();
  DiffBetweenSumAndWeights();
  SolveEnergyFunctional();
}

// Solve the energy functional for soft seeds
void SeededSegmentationSoft::SolveEnergyFunctional() {
  // Sparse matrix representing the linear system
  SparseMatrixType sparseA{m_sparseIs};
  sparseA += m_sparseL * m_sparseL;

  Eigen::CholmodDecomposition<SparseMatrixType, Eigen::Lower> sparseSolver;
  sparseSolver.setMode(Eigen::CholmodSupernodalLLt);
  sparseSolver.compute(sparseA);

  // Convert matrix of seeds into a vector
  auto const Matrix2Vector{[](Eigen::MatrixXd const &M) {
    auto const rows{M.rows()};
    auto const cols{M.cols()};
    Eigen::VectorXd V(rows * cols);
    for (auto i{0}; i < rows; ++i) {
      for (auto j{0}; j < cols; ++j) {
        V(i * cols + j) = M(i, j);
      }
    }
    return V;
  }};
  Eigen::VectorXd const seedsVector{Matrix2Vector(m_seeds)};

  // Solve
  Eigen::VectorXd const solutionVector{sparseSolver.solve(seedsVector)};

  // Convert the solution vector to matrix
  auto const Vector2Matrix{[](Eigen::VectorXd const &V, int rows, int cols) {
    Eigen::MatrixXd M(rows, cols);
    for (auto i{0}; i < rows; ++i) {
      for (auto j{0}; j < cols; ++j) {
        M(i, j) = V(i * cols + j);
      }
    }
    return M;
  }};
  m_solutionMatrix = Vector2Matrix(solutionVector, m_rows, m_cols);
}

// Create QImage using thresholded image as mask
void SeededSegmentationSoft::GetMaskedImage(cv::Mat const &mask, QImage &output,
                                            bool asBinary) const {
  auto const width{m_input->width()};
  auto const height{m_input->height()};

  for (auto i{0}; i < height; ++i) {
    for (auto j{0}; j < width; ++j) {
      QColor color;
      QPoint const qpt{j, i};
      if (mask.at<unsigned char>(cv::Point(j, i)) == 255)
        color = asBinary ? QColor(255, 255, 255, 255) : m_input->pixel(qpt);
      else
        color = asBinary ? QColor(0, 0, 0, 255) : QColor(0, 0, 0, 0);
      output.setPixelColor(qpt, color);
    }
  }
}

void SeededSegmentationSoft::GetOutput(QImage &output, bool asBinary) const {
  auto const width{m_input->width()};
  auto const height{m_input->height()};

  // Image from solution matrix
  cv::Mat grayscale(height, width, CV_8UC1);
  for (auto i{0}; i < m_rows; ++i) {
    for (auto j{0}; j < m_cols; ++j) {
      auto value{((m_solutionMatrix(i, j) + 1.0) / 2.0) * 255.0};

      value = std::clamp(value, 0.0, 255.0);
      if (std::isnan(value))
        value = 0.0;

      grayscale.at<unsigned char>(cv::Point(j, i)) =
          static_cast<unsigned char>(value);
    }
  }

  auto const mask{ApplyThresholding(grayscale)};

  GetMaskedImage(mask, output, asBinary);
}

} // namespace lc
