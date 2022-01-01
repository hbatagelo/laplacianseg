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

#include "seededsegmentationhard.h"

#include <Eigen/CholmodSupport>
#include <opencv2/opencv.hpp>

namespace lc {

SeededSegmentationHard::SeededSegmentationHard(QImage const *input,
                                               QImage const *seeds,
                                               QColor foreground,
                                               QColor background, double beta)
    : m_input{input},
      m_cols{input->width()},
      m_rows{input->height()},
      m_numPixels{m_rows * m_cols},
      m_beta{beta} {
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
  m_idxForeground.reserve(m_numPixels);
  m_idxBackground.reserve(m_numPixels);
  m_idxUnknown.reserve(m_numPixels);
  m_seeds = Eigen::MatrixXi(m_rows, m_cols);

  for (auto idx{0}, y{0}; y < m_rows; ++y) {
    for (auto x{0}; x < m_cols; ++x) {
      auto const color{seeds->pixelColor(x, y)};
      if (color == foreground) {
        m_seeds(y, x) = 1;
        m_idxBackground.push_back(idx);
      } else if (color == background) {
        m_seeds(y, x) = 0;
        m_idxForeground.push_back(idx);
      } else {
        m_seeds(y, x) = -1;
        m_idxUnknown.push_back(idx);
      }
      ++idx;
    }
  }
}

// Initialize the sparse matrices with the number of expected values in each row
// and set up the variables which will help us to compute the graph weight and
// the sum of the rows
void SeededSegmentationHard::SetUpGraphWeightAndSum() {
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

      if (!inBounds(pJ, boundsLo, boundsHi)) continue;

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

// Compute difference between the row sum of the weights and the weights
void SeededSegmentationHard::DiffBetweenSumAndWeights() {
  m_sparseL = m_sparseD;
  m_sparseL -= m_sparseWeights;
}

void SeededSegmentationHard::Compute() {
  SetUpGraphWeightAndSum();
  DiffBetweenSumAndWeights();
  SolveEnergyFunctional();
}

// Solve the energy functional for hard seeds
void SeededSegmentationHard::SolveEnergyFunctional() {
  std::vector<int> seededIdx;
  seededIdx.reserve(m_idxForeground.size() + m_idxBackground.size());
  seededIdx.insert(seededIdx.end(), m_idxBackground.begin(),
                   m_idxBackground.end());
  seededIdx.insert(seededIdx.end(), m_idxForeground.begin(),
                   m_idxForeground.end());

  std::vector<int> permutationVector;
  permutationVector.reserve(seededIdx.size() + m_idxUnknown.size());
  permutationVector.insert(permutationVector.end(), seededIdx.begin(),
                           seededIdx.end());
  permutationVector.insert(permutationVector.end(), m_idxUnknown.begin(),
                           m_idxUnknown.end());

  Eigen::VectorXi const permutation{
      Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(
          permutationVector.data(), long(permutationVector.size()))};
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> const perm(
      permutation);

  SparseMatrixType spLtwisted;
  spLtwisted = m_sparseL.twistedBy(perm.transpose());

  auto const bRows{std::ssize(seededIdx)};
  SparseMatrixType const sparseLs{spLtwisted.topLeftCorner(bRows, bRows)};

  auto const bCols{m_sparseL.cols() - bRows};
  SparseMatrixType const sparseB{spLtwisted.block(0, bRows, bRows, bCols)};
  SparseMatrixType const sparseBt{sparseB.transpose()};
  SparseMatrixType const sparseLu{spLtwisted.block(bRows, bRows, bCols, bCols)};

  Eigen::VectorXd xs{bRows};
  xs.setConstant(1);
  xs.segment(0, std::ssize(m_idxBackground)).setConstant(0);

  SparseMatrixType const sparseLHard{sparseBt * sparseB + sparseLu * sparseLu};
  Eigen::VectorXd const bHard{-(sparseLu * sparseBt + sparseBt * sparseLs) *
                              xs};

  Eigen::CholmodDecomposition<SparseMatrixType, Eigen::Lower> sparseSolver;
  sparseSolver.setMode(Eigen::CholmodSupernodalLLt);
  sparseSolver.compute(sparseLHard);
  m_solutionVector = sparseSolver.solve(bHard);
}

// Create QImage using thresholded image as mask
void SeededSegmentationHard::GetMaskedImage(cv::Mat const &mask, QImage &output,
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

void SeededSegmentationHard::GetOutput(QImage &output, bool asBinary) const {
  auto const width{m_input->width()};
  auto const height{m_input->height()};

  // Image from solution matrix
  cv::Mat grayscale(height, width, CV_8UC1);
  for (auto i{0}; i < m_rows; ++i) {
    for (auto j{0}; j < m_cols; ++j) {
      double value{m_seeds(i, j) == 1 ? 255.0 : 0.0};
      grayscale.at<unsigned char>(cv::Point(j, i)) =
          static_cast<unsigned char>(value);
    }
  }

  for (auto index{0}; index < m_solutionVector.size(); ++index) {
    auto value{(1.0 - m_solutionVector.coeff(index)) * 255.0};

    value = std::clamp(value, 0.0, 255.0);
    if (std::isnan(value)) value = 0.0;

    auto const offset{m_idxUnknown[static_cast<std::size_t>(index)]};
    auto const rowPi{offset / m_cols};
    auto const colPi{offset % m_cols};
    grayscale.at<unsigned char>(cv::Point(colPi, rowPi)) =
        static_cast<unsigned char>(value);
  }

  auto const mask{ApplyThresholding(grayscale)};

  GetMaskedImage(mask, output, asBinary);
}

}  // namespace lc
