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

#include "seededspsegmentation.h"

#include <Eigen/CholmodSupport>
#include <opencv2/core/core.hpp>

#include "util.h"

namespace lc {

SeededSPSegmentation::SeededSPSegmentation(QImage const *input,
                                           QImage const *seeds,
                                           int superpixelSize,
                                           double compactness)
    : m_input{input}, m_seeds{seeds}, m_slic{input, superpixelSize,
                                             compactness} {
  Q_ASSERT(input != nullptr);
  Q_ASSERT(seeds != nullptr);
}

void SeededSPSegmentation::Compute(bool useHardSeeds, QColor foreground,
                                   QColor background, double betaFeatures,
                                   double betaDistance) {
  m_useHardSeeds = useHardSeeds;

  m_slic.ComputeLabels();
  m_slic.ComputeGraph(m_seeds, foreground, background);
  m_slic.ComputeFeatures();

  CheckIfInSeeds();

  if (m_useHardSeeds) {
    std::vector<int> idxForeground;
    std::vector<int> idxBackground;
    std::vector<int> idxUnknown;
    idxForeground.reserve(m_slic.GetNumLabels());
    idxBackground.reserve(m_slic.GetNumLabels());
    idxUnknown.reserve(m_slic.GetNumLabels());

    // Get indices of seeded and non-seeded superpixels
    auto const &nodes{m_slic.GetLabelNodes()};
    for (auto i{0}; auto const &node : nodes) {
      switch (node.GetSeedType()) {
      case SeedType::Background:
        idxBackground.push_back(i);
        break;
      case SeedType::Foreground:
        idxForeground.push_back(i);
        break;
      case SeedType::None:
        idxUnknown.push_back(i);
        break;
      }
      ++i;
    }

    auto numColumns{std::ssize(nodes) -
                    (std::ssize(idxBackground) + std::ssize(idxForeground))};
    m_solutionVector = Eigen::VectorXd::Zero(numColumns);
  } else {
    m_solutionVector = Eigen::VectorXd::Zero(m_slic.GetNumLabels());
  }

  SetUpGraphWeightAndSum(betaFeatures, betaDistance);
  DiffBetweenSumAndWeights();

  if (m_useHardSeeds)
    SolveEnergyFunctionalHard();
  else
    SolveEnergyFunctionalSoft();
}

// Returns the beta which results in the largest variation of weights.
// From betaMin to betaMax (using betaStep as step), returns the beta which
// results in the greatest variation of weights.
// If global is false, the greatest variation of weights is computed w.r.t.
// the neighborhood of each superpixel, using a histogram of occurrences of
// each beta.
double SeededSPSegmentation::ComputeDistanceBeta(double betaMin, double betaMax,
                                                 double betaStep, bool global) {
  auto const &nodes{m_slic.GetLabelNodes()};
  auto beta{0.0};
  auto bestBeta{0.0};

  // Compute minimum and maximum neighbors weights (minWeight, maxWeight) of
  // a given superpixel (node).
  auto ComputeWeightsFromNeighbors{
      [&nodes, &beta](auto const &node, auto &minWeight, auto &maxWeight) {
        auto const &centerOfMass{node.GetCenterOfMass()};
        for (auto const &neighbor : node.GetNeighbors()) {
          auto const &centerOfMassNeighbor{nodes[neighbor].GetCenterOfMass()};
          // Distance between centroids
          auto const distance{
              Eigen::Vector2d{centerOfMassNeighbor - centerOfMass}.norm()};
          auto const weight{std::exp(-beta * distance)};
          minWeight = std::min(minWeight, weight);
          maxWeight = std::max(maxWeight, weight);
        }
      }};

  // Global beta
  if (global) {
    auto maxGlobalWeightDiff{-1.0};

    beta = betaMin;
    while (beta < betaMax) {
      auto minWeight{std::numeric_limits<double>::max()};
      auto maxWeight{std::numeric_limits<double>::lowest()};

      // For each superpixel node
      for (auto const &node : nodes) {
        ComputeWeightsFromNeighbors(node, minWeight, maxWeight);
      }

      // Keep global maximum
      auto const curWeightDiff{maxWeight - minWeight};
      if (curWeightDiff > maxGlobalWeightDiff) {
        maxGlobalWeightDiff = curWeightDiff;
        bestBeta = beta;
      }
      beta += betaStep;
    }
  }
  // Local (neighborhood) beta
  else {
    auto numHistogramBins{
        static_cast<std::size_t>(std::ceil((betaMax - betaMin) / betaStep))};
    std::vector<int> localBetaHistogram(numHistogramBins, 0);

    // For each superpixel node
    for (auto const &node : nodes) {
      auto maxLocalWeightDiff{-1.0};
      auto bestLocalBeta{0.0};

      beta = betaMin;
      while (beta < betaMax) {
        auto minWeight{std::numeric_limits<double>::max()};
        auto maxWeight{std::numeric_limits<double>::lowest()};

        ComputeWeightsFromNeighbors(node, minWeight, maxWeight);

        // Keep local maximum
        auto const curWeightDiff{maxWeight - minWeight};
        if (curWeightDiff > maxLocalWeightDiff) {
          maxLocalWeightDiff = curWeightDiff;
          bestLocalBeta = beta;
        }
        beta += betaStep;
      }

      auto const histogramIndex{
          static_cast<std::size_t>((bestLocalBeta - betaMin) / betaStep)};
      localBetaHistogram[histogramIndex]++;
    }

    auto const maxHistogramIndex{
        std::distance(localBetaHistogram.begin(),
                      std::ranges::max_element(localBetaHistogram))};
    bestBeta = static_cast<double>(maxHistogramIndex) * betaStep + betaMin;
  }

  return bestBeta;
}

// Initialize the sparse matrices with the number of expected values in each
// row and instantiate the variables which will help us to compute the graph
// weight and the sum of the rows
void SeededSPSegmentation::SetUpGraphWeightAndSum(double betaFeatures,
                                                  double betaDistance) {
  constexpr auto epsilon{1e-6};
  constexpr auto epsilonInv{1.0 / epsilon};

  auto const n{m_slic.GetNumLabels()};
  auto const maxValency{m_slic.GetMaxVertexValency()};

  m_sparseWeights = Eigen::SparseMatrix<double>(n, n);
  m_sparseWeights.reserve(Eigen::VectorXi::Constant(n, maxValency));

  m_sparseD = Eigen::SparseMatrix<double>(n, n);
  m_sparseD.reserve(Eigen::VectorXi::Constant(n, 1));

  if (betaDistance < 0)
    betaDistance = ComputeDistanceBeta();

  // For each superpixel node
  auto const &nodes{m_slic.GetLabelNodes()};
  for (auto idx{0}; auto const &node : nodes) {
    auto summ{0.0};

    auto const &featureVector{node.GetFeatureVector()};
    auto const &centroidCurrent{node.GetCenterOfMass()};

    // For each neighbor superpixel
    for (auto const &neighIdx : node.GetNeighbors()) {
      auto wij{m_sparseWeights.coeff(idx, neighIdx)};

      // Still untouched?
      if (almostEquals(wij, 0.0)) {
        // Euclidean distance between centroids
        auto const &centroidNeighbor{nodes[neighIdx].GetCenterOfMass()};
        auto const distance{
            Eigen::Vector2d{centroidNeighbor - centroidCurrent}.norm()};

        // Distance from Lab attributes
        auto const &featureVectorNeighbor{nodes[neighIdx].GetFeatureVector()};
        Eigen::VectorXd labDiff{featureVector - featureVectorNeighbor};

        // Weight
        auto const weight{labDiff.size() > 0 ? labDiff = labDiff.cwiseAbs(),
                          labDiff.lpNorm<Eigen::Infinity>() : 0.0};
        Q_ASSERT(!std::isnan(weight));
        wij = std::exp(-betaFeatures * weight) *
              std::exp(-betaDistance * distance);
        wij = std::round(wij * epsilonInv) * epsilon;
        wij += epsilon;

        m_sparseWeights.insert(idx, neighIdx) = wij;
        m_sparseWeights.insert(neighIdx, idx) = wij;
      }

      summ += wij;
    }

    // Insert the sum in the diagonal matrix
    m_sparseD.insert(idx, idx) = summ;
    ++idx;
  }

  m_sparseWeights.makeCompressed();
  m_sparseD.makeCompressed();
}

// Set up the diagonal sparse matrix with ones if the pixel belongs to the
// background or the foreground, and zeros otherwise
void SeededSPSegmentation::CheckIfInSeeds() {
  if (m_useHardSeeds)
    return;

  auto const n{m_slic.GetNumLabels()};
  m_sparseIs = Eigen::SparseMatrix<double>(n, n);
  m_sparseIs.reserve(Eigen::VectorXi::Constant(n, 1));
  m_seedsVector = Eigen::VectorXd(n);

  for (auto idx{0}; auto const &node : m_slic.GetLabelNodes()) {
    auto type{node.GetSeedType()};
    m_seedsVector(idx) = static_cast<int>(type);

    if (type != SeedType::None)
      m_sparseIs.insert(idx, idx) = 1;

    ++idx;
  }

  m_sparseIs.makeCompressed();
}

// Compute difference between the row sum of the weights and the weights
void SeededSPSegmentation::DiffBetweenSumAndWeights() {
  m_sparseL = m_sparseD;
  m_sparseL -= m_sparseWeights;
}

// Rewrite the energy functional for soft seeds
void SeededSPSegmentation::SolveEnergyFunctionalSoft() {
  SparseMatrixType sparseA{m_sparseIs};
  sparseA += m_sparseL * m_sparseL;

  Eigen::CholmodDecomposition<SparseMatrixType, Eigen::Lower> sparseSolver;
  sparseSolver.setMode(Eigen::CholmodSupernodalLLt);
  sparseSolver.compute(sparseA);

  Eigen::VectorXd solutionVector{sparseSolver.solve(m_seedsVector)};
  m_solutionVector += solutionVector;
}

// Solve the energy functional for hard seeds
void SeededSPSegmentation::SolveEnergyFunctionalHard() {
  auto const n{m_slic.GetNumLabels()};

  // Get indices of seeded and non-seeded superpixels
  std::vector<int> idxForeground;
  std::vector<int> idxBackground;
  std::vector<int> idxUnknown;
  idxForeground.reserve(n);
  idxBackground.reserve(n);
  idxUnknown.reserve(n);
  for (auto idx{0}; auto const &node : m_slic.GetLabelNodes()) {
    auto type{node.GetSeedType()};
    switch (type) {
    case SeedType::Background:
      idxBackground.push_back(idx);
      break;
    case SeedType::Foreground:
      idxForeground.push_back(idx);
      break;
    case SeedType::None:
      idxUnknown.push_back(idx);
      break;
    }
    ++idx;
  }

  if (idxUnknown.empty())
    return; // Nothing to solve

  std::vector<int> seededIdx;
  seededIdx.reserve(idxForeground.size() + idxBackground.size());
  seededIdx.insert(seededIdx.end(), idxBackground.begin(), idxBackground.end());
  seededIdx.insert(seededIdx.end(), idxForeground.begin(), idxForeground.end());

  std::vector<int> permutationVector;
  permutationVector.reserve(seededIdx.size() + idxUnknown.size());
  permutationVector.insert(permutationVector.end(), seededIdx.begin(),
                           seededIdx.end());
  permutationVector.insert(permutationVector.end(), idxUnknown.begin(),
                           idxUnknown.end());

  Eigen::VectorXi const permutation{
      Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(
          permutationVector.data(), std::ssize(permutationVector))};
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> const perm{
      permutation};

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
  xs.segment(0, std::ssize(idxBackground)).setConstant(0);

  SparseMatrixType const sparseLHard{sparseBt * sparseB + sparseLu * sparseLu};
  Eigen::VectorXd const bHard{-(sparseLu * sparseBt + sparseBt * sparseLs) *
                              xs};

  Eigen::CholmodDecomposition<SparseMatrixType, Eigen::Lower> sparseSolver;
  sparseSolver.setMode(Eigen::CholmodSupernodalLLt);
  sparseSolver.compute(sparseLHard);
  m_solutionVector = sparseSolver.solve(bHard);
}

namespace {
// Return true if pixel was found in image; false otherwise
// segImage is the segmented image. This is needed for soft seeds to guarantee
// that the flood fill seed is inside a foreground area.
bool GetFirstNonZeroPixel(cv::Mat const &image, cv::Mat const &segImage,
                          cv::Point &point, bool useHardSeeds) {
  for (auto i{0}; i < image.rows; ++i) {
    for (auto j{0}; j < image.cols; ++j) {
      cv::Point const cpt{j, i};
      if (image.at<unsigned char>(cpt) != 0) {
        if (!useHardSeeds) {
          if (segImage.at<unsigned char>(cpt) == 0)
            continue;
        }
        point = cpt;
        return true;
      }
    }
  }
  return false;
}
} // namespace

void SeededSPSegmentation::RemoveUnlabeledConnectedComponents(
    cv::Mat const &seeds, cv::Mat &mask) const {
  // Flood fill to remove components without seeds
  cv::Mat imUnusedSeeds;
  seeds.copyTo(imUnusedSeeds);

  // Repeat until there are no seeds left
  cv::Point seedPoint;
  while (GetFirstNonZeroPixel(imUnusedSeeds, mask, seedPoint, m_useHardSeeds)) {
    cv::floodFill(mask, seedPoint, 128);

    cv::Mat imFilled;
    cv::threshold(mask, imFilled, 128, 0, cv::THRESH_TOZERO_INV);
    cv::threshold(imFilled, imFilled, 0, 255, cv::THRESH_BINARY);

    imUnusedSeeds = imUnusedSeeds - imFilled;
  }

  cv::threshold(mask, mask, 128, 0, cv::THRESH_TOZERO_INV); // Isolate 128
  cv::threshold(mask, mask, 0, 255, cv::THRESH_BINARY);     // Map to 255
}

// Create QImage using thresholded image as mask
void SeededSPSegmentation::GetMaskedImage(cv::Mat const &mask, QImage &output,
                                          bool asBinary) const {
  auto const width{m_input->width()};
  auto const height{m_input->height()};

  for (auto i{0}; i < height; ++i) {
    for (auto j{0}; j < width; ++j) {
      QColor color;
      QPoint const qpt{j, i};

      if (mask.at<unsigned char>(cv::Point(j, i)) == 255)
        color = asBinary ? QColor(255, 255, 255) : m_input->pixel(qpt);
      else
        color = asBinary ? QColor(0, 0, 0, 255) : QColor(0, 0, 0, 0);

      output.setPixelColor(qpt, color);
    }
  }
}

void SeededSPSegmentation::GetOutput(QImage &output, bool asBinary) const {
  auto const width{m_input->width()};
  auto const height{m_input->height()};

  // Image from solution vector
  cv::Mat grayscale(height, width, CV_8UC1);
  for (auto solutionIdx{0}; auto const &node : m_slic.GetLabelNodes()) {
    auto value{0.0};

    if (m_useHardSeeds) {
      switch (node.GetSeedType()) {
      case SeedType::Background:
        value = 0.0;
        break;
      case SeedType::Foreground:
        value = 1.0;
        break;
      case SeedType::None:
        value = m_solutionVector(solutionIdx);
        solutionIdx++;
        break;
      }
      value *= 255.0;
    } else {
      value = ((m_solutionVector(solutionIdx) + 1.0) / 2.0) * 255;
      solutionIdx++;
    }

    value = std::clamp(value, 0.0, 255.0);
    if (std::isnan(value))
      value = 0.0;

    // Pixel coordinates for this superpixel
    for (auto const &pt : node.GetInnerPixels()) {
      grayscale.at<unsigned char>(cv::Point(pt.x(), pt.y())) =
          static_cast<unsigned char>(value);
    }
  }

  // Foreground and background seeds
  cv::Mat seedsFg{cv::Mat::zeros(grayscale.size(), grayscale.type())};
  cv::Mat seedsBg{cv::Mat::zeros(grayscale.size(), grayscale.type())};

  // Copy seeds to imSeedsFg and imSeedsBg
  for (auto const &node : m_slic.GetLabelNodes()) {
    auto const type{node.GetSeedType()};
    if (type == SeedType::None)
      continue;

    auto &seeds{type == SeedType::Foreground ? seedsFg : seedsBg};
    for (auto const &pt : node.GetInnerPixels()) {
      seeds.at<unsigned char>(cv::Point(pt.x(), pt.y())) = 255;
    }
  }

  auto mask{ApplyThresholding(grayscale)};

  RemoveUnlabeledConnectedComponents(seedsFg, mask);

  GetMaskedImage(mask, output, asBinary);
}

} // namespace lc
