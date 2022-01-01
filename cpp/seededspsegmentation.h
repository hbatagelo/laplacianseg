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

#ifndef SEEDEDSPSEGMENTATION_H
#define SEEDEDSPSEGMENTATION_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <QImage>
#include <stdexcept>

#include "opencv2/opencv.hpp"
#include "slic.h"

namespace lc {

class SeededSPSegmentation {
  using SparseMatrixType = Eigen::SparseMatrix<double>;

  QImage const *m_input;  // Input image
  QImage const *m_seeds;  // Input seeds image

  SLIC m_slic;

  // Seeds vector (single value per superpixel)
  Eigen::VectorXd m_seedsVector;

  // Sparse matrix storing the weights
  SparseMatrixType m_sparseWeights;

  // Sparse matrices minimizing the energy functional
  SparseMatrixType m_sparseD;
  SparseMatrixType m_sparseIs;
  SparseMatrixType m_sparseL;

  Eigen::VectorXd m_solutionVector;

  bool m_useHardSeeds{};

  double ComputeDistanceBeta(double betaMin = 0.1, double betaMax = 0.5,
                             double betaStep = 0.01, bool global = false);
  void SetUpGraphWeightAndSum(double betaFeatures, double betaDistance);
  void CheckIfInSeeds();
  void DiffBetweenSumAndWeights();
  void SolveEnergyFunctionalSoft();
  void SolveEnergyFunctionalHard();
  void RemoveUnlabeledConnectedComponents(cv::Mat const &seeds,
                                          cv::Mat &mask) const;
  void GetMaskedImage(cv::Mat const &mask, QImage &output, bool asBinary) const;

 public:
  SeededSPSegmentation(QImage const *input, QImage const *seeds,
                       int superpixelSize, double compactness);
  void Compute(bool useHardSeeds = false, QColor foreground = Qt::red,
               QColor background = Qt::blue, double betaFeatures = 88.0,
               double betaDistance = -1.0);
  void GetOutput(QImage &output, bool asBinary = false) const;
};

}  // namespace lc

#endif  // SEEDEDSPSEGMENTATION_H
