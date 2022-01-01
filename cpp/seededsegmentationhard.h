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

#ifndef SEEDEDSEGMENTATIONHARD_H
#define SEEDEDSEGMENTATIONHARD_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <QImage>
#include <opencv2/opencv.hpp>

#include "util.h"

namespace lc {

class SeededSegmentationHard {
  using SparseMatrixType = Eigen::SparseMatrix<double>;

  QImage const *m_input;
  int const m_cols;
  int const m_rows;
  int const m_numPixels;
  double const m_beta;

  // Input image matrix (3 channels)
  std::array<Eigen::MatrixXd, 3> m_channels{};

  // Seeds matrix (single value per pixel)
  Eigen::MatrixXi m_seeds;

  // Linear indices of foreground, background and unknown seeds
  std::vector<int> m_idxForeground, m_idxBackground, m_idxUnknown;

  // Sparse matrix storing the weights
  SparseMatrixType m_sparseWeights;

  // Sparse matrices minimizing the energy functional
  SparseMatrixType m_sparseD;
  SparseMatrixType m_sparseIs;
  SparseMatrixType m_sparseL;

  Eigen::VectorXd m_solutionVector;

  void SetUpGraphWeightAndSum();
  void DiffBetweenSumAndWeights();
  void SolveEnergyFunctional();
  void GetMaskedImage(cv::Mat const &mask, QImage &output, bool asBinary) const;

public:
  SeededSegmentationHard(QImage const *input, QImage const *seeds,
                         QColor foreground = Qt::red,
                         QColor background = Qt::blue, double beta = 500.0);
  void Compute();
  void GetOutput(QImage &output, bool asBinary = false) const;
};

} // namespace lc

#endif // SEEDEDSEGMENTATIONHARD_H
