/****************************************************************************
**
** Copyright (C) 2020 Harlen Batagelo <hbatagelo@gmail.com> and
**                    João Paulo Gois <jpgois@gmail.com>.
**
** This file is part of the implementation of the paper
** 'Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation'
** by Wallace Casaca, João Paulo Gois, Harlen Batagelo, Gabriel Taubin and
** Luis Gustavo Nonato.
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

#include <QVector2D>

#include "opencv2/core/core.hpp"
#include "seededspsegmentation.h"
#include "util.h"

// Lp-norm used with the feature vector
const int featureNorm = Eigen::Infinity;

SeededSPSegmentation::SeededSPSegmentation(QImage *inputImage,
                                           QImage *seedsImage,
                                           int superpixelSize,
                                           double compactness)
    : m_image(inputImage),
      m_seeds(seedsImage),
      m_slic(inputImage, superpixelSize, compactness)
{
}

bool SeededSPSegmentation::Compute(bool useHardSeeds,
                                   QColor fgColor,
                                   QColor bgColor,
                                   double beta,
                                   double betaDistance)
{
    m_useHardSeeds = useHardSeeds;

    m_slic.ComputeLabels();
    // m_seeds is used to compute the number of seeded pixels inside each superpixel
    m_slic.ComputeGraph(m_seeds, fgColor, bgColor);
    m_slic.ComputeFeatures();

    CheckIfInSeedsOrNot();

    if (m_useHardSeeds)
    {
        // Get indices of seeded and non-seeded superpixels
        std::vector<int> idxForeground, idxBackground, idxUnknown;
        auto &nodes = m_slic.GetLabelNodes();
        for (auto i = 0; i < nodes.size(); ++i)
        {
            int seedType = nodes[i].GetSeedType();
            if (seedType == -1) // Background
                idxBackground.push_back(i);
            else if (seedType == 1) // Foreground
                idxForeground.push_back(i);
            else // Not seeded
                idxUnknown.push_back(i);
        }

        std::vector<int> seededidx;
        seededidx.insert(seededidx.end(), idxBackground.begin(), idxBackground.end());
        seededidx.insert(seededidx.end(), idxForeground.begin(), idxForeground.end());

        auto b_cols = nodes.size() - seededidx.size();
        m_solutionVector = Eigen::VectorXd::Zero(b_cols);
    }
    else
    {
        m_solutionVector = Eigen::VectorXd::Zero(m_slic.GetNumLabels());
    }

    SetUpGraphWeightAndSum(beta, betaDistance);
    DiffBetweenSumAndWeights();

    if (m_useHardSeeds)
        SolveHardLinearSystem();
    else
        SolveSoftLinearSystem();

    return true;
}

void SeededSPSegmentation::GetOutput(QImage &outputImage, bool asBinary)
{
    auto &nodes = m_slic.GetLabelNodes();

    int width = m_image->width();
    int height = m_image->height();

    bool otsu = true;
    // If otsuVariant is true, the threshold is
    // given by Isol > (1-otsu) instead of Isol > otsu
    bool otsuVariant = true;

    bool removeConnectedComponentsWithoutLabels = true;

    cv::Mat imBw(height, width, CV_8UC1);
    cv::Mat imGray(height, width, CV_8UC1);
    cv::Mat imSeedsFg(height, width, CV_8UC1); // Foreground seeds
    cv::Mat imSeedsBg(height, width, CV_8UC1); // Background seeds

    // Mount solution vector as an OpenCV matrix
    auto solutionIdx = 0;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        double value = 0;

        if (m_useHardSeeds)
        {
            int seedType = nodes[i].GetSeedType();
            if (seedType == 1) // Foreground
            {
                value = 1;
            }
            else if (seedType == -1) // Background
            {
                value = 0;
            }
            else
            {
                value = m_solutionVector(solutionIdx);
                solutionIdx++;
            }
            value *= 255.0;
        }
        else
        {
            value = m_solutionVector(solutionIdx);
            value = ((value + 1) / 2.0) * 255;
            solutionIdx++;
        }

        value = std::clamp(value, 0.0, 255.0);
        if (std::isnan(value))
        {
            value = 0;
        }

        // List of (x, y) pixel coordinates for this superpixel
        for (const auto &pt : nodes[i].GetPixels())
        {
            auto cpt = cv::Point(pt.x(), pt.y());
            imGray.at<uchar>(cpt) = uchar(value);
        }
    }

    imSeedsFg.setTo(0);
    imSeedsBg.setTo(0);

    // Copy seeds to imSeedsFg and imSeedsBg
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        int seedType = nodes[i].GetSeedType();
        if (seedType == 1) // Foreground
        {
            for (const auto &pt : nodes[i].GetPixels())
            {
                auto cpt = cv::Point(pt.x(), pt.y());
                imSeedsFg.at<uchar>(cpt) = 255;
            }
        }
        else if (seedType == -1) // Background
        {
            for (const auto &pt : nodes[i].GetPixels())
            {
                auto cpt = cv::Point(pt.x(), pt.y());
                imSeedsBg.at<uchar>(cpt) = 255;
            }
        }
    }

    // Threshold with Otsu
    if (otsu)
    {
        if (otsuVariant)
        {
            double otsuThreshold = cv::threshold(imGray, imBw, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
            if (almostEquals(otsuThreshold, 0.0))
                otsuThreshold = 255;
            cv::threshold(imGray, imBw, 255 - otsuThreshold, 255, CV_THRESH_BINARY);
        }
        else
        {
            cv::threshold(imGray, imBw, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
        }
    }
    else
    {
        cv::threshold(imGray, imBw, 127, 255, CV_THRESH_BINARY);
    }

    if (removeConnectedComponentsWithoutLabels)
    {
        // Flood fill to remove components without seeds
        cv::Mat imUnusedSeeds;
        imSeedsFg.copyTo(imUnusedSeeds);

        // Repeat this until there are no seeds left
        cv::Point seedPoint;
        while (GetFirstNonZeroPixel(imUnusedSeeds, imBw, seedPoint))
        {
            cv::floodFill(imBw, seedPoint, 128);

            cv::Mat imFilled;
            cv::threshold(imBw, imFilled, 128, 0, cv::THRESH_TOZERO_INV); // Isolate 128
            cv::threshold(imFilled, imFilled, 0, 255, cv::THRESH_BINARY); // Map to 255

            imUnusedSeeds = imUnusedSeeds - imFilled;
        }

        cv::threshold(imBw, imBw, 128, 0, cv::THRESH_TOZERO_INV); // Isolate 128
        cv::threshold(imBw, imBw, 0, 255, cv::THRESH_BINARY); // Map to 255
    }

    // Copy thresholded image to output QImage
    for (auto i = 0; i < height; ++i)
    {
        for (auto j = 0; j < width; ++j)
        {
            auto cpt = cv::Point(j, i);
            auto pt = QPoint(j, i);
            if (imBw.at<uchar>(cpt) == 255)
            {
                if (asBinary)
                    outputImage.setPixelColor(pt, QColor(255,255,255));
                else
                    outputImage.setPixelColor(pt, m_image->pixel(pt));
            }
            else
            {
                if (asBinary)
                    outputImage.setPixelColor(pt, QColor(0,0,0,255));
                else
                    outputImage.setPixelColor(pt, QColor(0,0,0,0));
            }
        }
    }
}

// Return true if pixel was found in image; false otherwise
// segImage is the segmented image. This is needed for soft seeds, to guarantee
// that the flood fill seed is inside a foreground area.
bool SeededSPSegmentation::GetFirstNonZeroPixel(cv::Mat &image,
                                                cv::Mat &segImage,
                                                cv::Point &point)
{
    for (auto i = 0; i < image.rows; ++i)
    {
        for (auto j = 0; j < image.cols; ++j)
        {
            cv::Point cpt = cv::Point(j, i);

            if (image.at<uchar>(cpt) != 0)
            {
                if (!m_useHardSeeds)
                {
                    if (segImage.at<uchar>(cpt) == 0)
                        continue;
                }
                point = cpt;
                return true;
            }
        }
    }
    return false;
}

// From betaMin to betaMax (using betaStep as step), returns the beta which
// results in the greatest variation of weights.
// If global is false, the greatest variation of weights is computed w.r.t.
// the neighborhood of each superpixel, using a histogram of occurrences of
// each beta.
double SeededSPSegmentation::ComputeDistanceBeta(double betaMin,
                                                 double betaMax,
                                                 double betaStep,
                                                 bool global)
{
    auto maxGlobalWeightDiff = -1.0;
    auto bestBeta = 0.0;

    auto &nodes = m_slic.GetLabelNodes();

    // Global beta
    if (global)
    {
        for (double beta = betaMin; beta < betaMax; beta += betaStep)
        {
            auto minGlobalWeight = std::numeric_limits<double>::max();
            auto maxGlobalWeight = std::numeric_limits<double>::lowest();

            // For each superpixel node
            for (auto i = 0; i < nodes.size(); ++i)
            {
                const auto& centroidCurrent = nodes[i].GetCenterOfMass();

                // For each neighboring superpixel node
                const auto &neighbors = nodes[i].GetNeighbors();
                for (auto j = 0; j < neighbors.size(); ++j)
                {
                    // Compute Euclidean distance between centroids
                    const auto& centroidNeighbor = nodes[neighbors[j]].GetCenterOfMass();
                    auto distance = QVector2D(centroidNeighbor - centroidCurrent).length();

                    auto weight = std::exp(-beta * distance);
                    minGlobalWeight = std::min(minGlobalWeight, weight);
                    maxGlobalWeight = std::max(maxGlobalWeight, weight);
                }
            }

            // Keep global maximum
            auto curWeightDiff = maxGlobalWeight - minGlobalWeight;
            if (curWeightDiff > maxGlobalWeightDiff)
            {
                maxGlobalWeightDiff = curWeightDiff;
                bestBeta = beta;
            }
        }
    }
    // Local (neighborhood) beta
    else
    {
        size_t numHistogramBins = std::ceil((betaMax - betaMin) / betaStep);
        std::vector<int> localBetaHistogram;
        localBetaHistogram.resize(numHistogramBins, 0);

        // For each superpixel node
        for (auto i = 0; i < nodes.size(); ++i)
        {
            auto maxLocalWeightDiff = -1.0;
            auto bestLocalBeta = 0.0;

            const auto& centroidCurrent = nodes[i].GetCenterOfMass();

            for (auto beta = betaMin; beta < betaMax; beta += betaStep)
            {
                auto minLocalWeight = std::numeric_limits<double>::max();
                auto maxLocalWeight = std::numeric_limits<double>::lowest();

                // For each neighboring superpixel node
                const auto &neighbors = nodes[i].GetNeighbors();
                for (auto j = 0; j < neighbors.size(); ++j)
                {
                    // Compute Euclidean distance between centroids
                    const auto& centroidNeighbor = nodes[neighbors[j]].GetCenterOfMass();
                    auto distance = QVector2D(centroidNeighbor - centroidCurrent).length();

                    auto weight = std::exp(-beta * distance);
                    minLocalWeight = std::min(minLocalWeight, weight);
                    maxLocalWeight = std::max(maxLocalWeight, weight);
                }

                // Keep local maximum
                double curWeightDiff = maxLocalWeight - minLocalWeight;
                if (curWeightDiff > maxLocalWeightDiff)
                {
                    maxLocalWeightDiff = curWeightDiff;
                    bestLocalBeta = beta;
                }
            }

            size_t histogramIndex = (bestLocalBeta - betaMin) / betaStep;
            Q_ASSERT(histogramIndex < numHistogramBins);
            localBetaHistogram[histogramIndex]++;
        }

        // Determine which bin has the most votes
        auto maxHistogramValue = 0;
        auto maxHistogramIndex = 0;
        for (auto i = 0; i < numHistogramBins; ++i)
        {
            if (localBetaHistogram[i] > maxHistogramValue)
            {
                maxHistogramValue = localBetaHistogram[i];
                maxHistogramIndex = i;
            }
        }

        bestBeta = maxHistogramIndex * betaStep + betaMin;
    }

    return bestBeta;
}

// Initialize the Sparse Matrices with the number of expected values in
// each row and instantiate the variables which will help us to compute
// the graph weight and the sum of the rows
void SeededSPSegmentation::SetUpGraphWeightAndSum(double betaFeatures, double betaDistance)
{
    const double epsilon = 1e-6;

    m_sparseWeights = Eigen::SparseMatrix<double>(m_slic.GetNumLabels(), m_slic.GetNumLabels());
    m_sparseWeights.reserve(Eigen::VectorXi::Constant(m_slic.GetNumLabels(), m_slic.GetMaxVertexValency()));
    m_sparseD = Eigen::SparseMatrix<double>(m_slic.GetNumLabels(), m_slic.GetNumLabels());
    m_sparseD.reserve(Eigen::VectorXi::Constant(m_slic.GetNumLabels(), 1));

    auto &nodes = m_slic.GetLabelNodes();

    if (betaDistance < 0)
        betaDistance = ComputeDistanceBeta();

    // For each superpixel node
    for (auto i = 0; i < nodes.size(); ++i)
    {
        double summ = 0;

        const auto& featureVector = nodes[i].GetFeatureVector();
        const auto& centroidCurrent = nodes[i].GetCenterOfMass();

        // For each neighboring superpixel node
        const auto &neighbors = nodes[i].GetNeighbors();
        for (auto j = 0; j < neighbors.size(); ++j)
        {
            double wij = m_sparseWeights.coeff(i, neighbors[j]);

            if (almostEquals(wij, 0.0)) // Still untouched?
            {
                auto weight = 0.0;

                // Compute Euclidean distance between centroids
                const auto& centroidNeighbor = nodes[neighbors[j]].GetCenterOfMass();
                auto distance = QVector2D(centroidNeighbor - centroidCurrent).length();

                // Compute distance from Lab attributes
                const auto& featureVectorNeighbor = nodes[neighbors[j]].GetFeatureVector();
                Eigen::VectorXd labDiff = featureVector - featureVectorNeighbor;

                auto labDiffNorm = 0.0;
                if (labDiff.size() > 0)
                {
                    labDiff = labDiff.cwiseAbs();
                    labDiffNorm = labDiff.lpNorm<featureNorm>();
                }

                weight = labDiffNorm;

                Q_ASSERT(!std::isnan(weight));

                // Compute the weight (features)
                wij = std::exp(-betaFeatures * weight) * std::exp(-betaDistance * distance);
                wij = std::round(wij * 1e6) * 1e-6;
                wij += epsilon;

                m_sparseWeights.insert(i, neighbors[j]) = wij;
                m_sparseWeights.insert(neighbors[j], i) = wij;
            }

            summ += wij;
        }

        // Insert the sum in the diagonal matrix
        m_sparseD.insert(i, i) = summ;
    }
}

// Check whether each superpixel is a seed. If so, insert 1 in the diagonal of
// the sparse matrix; do nothing otherwise
void SeededSPSegmentation::CheckIfInSeedsOrNot()
{
    if (m_useHardSeeds)
        return;

    m_sparseIs = Eigen::SparseMatrix<double>(m_slic.GetNumLabels(), m_slic.GetNumLabels());
    m_sparseIs.reserve(Eigen::VectorXi::Constant(m_slic.GetNumLabels(), 1));
    m_seedsVector = Eigen::VectorXd(m_slic.GetNumLabels());

    auto &nodes = m_slic.GetLabelNodes();
    for (auto i = 0; i < nodes.size(); ++i)
    {
        int seedType = nodes[i].GetSeedType();
        m_seedsVector(i) = seedType;

        if (seedType != 0)
            m_sparseIs.insert(i, i) = 1;
    }
}

// Compute difference between the row sum of the weights and the weights
void SeededSPSegmentation::DiffBetweenSumAndWeights()
{
    m_sparseL = Eigen::SparseMatrix<double>(m_slic.GetNumLabels(), m_slic.GetNumLabels());
    m_sparseL.reserve(Eigen::VectorXi::Constant(m_slic.GetNumLabels(), m_slic.GetMaxVertexValency()));
    m_sparseL = m_sparseD - m_sparseWeights;
}

// Rewrite the energy function in matrices and minimize this one
// using Cholesky factorization
void SeededSPSegmentation::SolveSoftLinearSystem()
{
    typedef Eigen::SparseMatrix<double> SparseMatrixType;

    SparseMatrixType sparseA;

    sparseA = SparseMatrixType(m_slic.GetNumLabels(), m_slic.GetNumLabels());
    sparseA.reserve(Eigen::VectorXi::Constant(m_slic.GetNumLabels(), m_slic.GetMaxVertexValency()));

    sparseA = (m_sparseIs + m_sparseL * m_sparseL);

    Eigen::CholmodDecomposition<SparseMatrixType, Eigen::Lower> sparseSolver;
    sparseSolver.setMode(Eigen::CholmodSupernodalLLt);
    sparseSolver.compute(sparseA);

    Eigen::VectorXd solutionVector = sparseSolver.solve(m_seedsVector);
    m_solutionVector += solutionVector;
}

void SeededSPSegmentation::SolveHardLinearSystem()
{
    typedef Eigen::SparseMatrix<double> SparseMatrixType;

    // Get indices of seeded and non-seeded superpixels
    std::vector<int> idxForeground, idxBackground, idxUnknown;
    auto &nodes = m_slic.GetLabelNodes();
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        int seedType = nodes[i].GetSeedType();
        if (seedType == -1) // Background
            idxBackground.push_back(i);
        else if (seedType == 1) // Foreground
            idxForeground.push_back(i);
        else // Not seeded
            idxUnknown.push_back(i);
    }

    if (idxUnknown.size() == 0)
        return; // Nothing to solve

    std::vector<int> seededidx;
    seededidx.insert(seededidx.end(), idxBackground.begin(), idxBackground.end());
    seededidx.insert(seededidx.end(), idxForeground.begin(), idxForeground.end());

    std::vector<int> vectorpermutation;
    vectorpermutation.insert(vectorpermutation.end(), seededidx.begin(), seededidx.end());
    vectorpermutation.insert(vectorpermutation.end(), idxUnknown.begin(), idxUnknown.end());

    Eigen::VectorXi permutation;
    permutation = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(vectorpermutation.data(), vectorpermutation.size());
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,int> perm(permutation);

    SparseMatrixType spLtwisted;
    spLtwisted = m_sparseL.twistedBy(perm.transpose());

    auto b_rows = seededidx.size();

    SparseMatrixType sparseLs;
    sparseLs = spLtwisted.topLeftCorner(b_rows, b_rows);

    auto b_cols = m_sparseL.cols() - b_rows;
    SparseMatrixType sparseB;
    sparseB= spLtwisted.block(0, b_rows, b_rows, b_cols);
    SparseMatrixType sparseBt = sparseB.transpose();

    SparseMatrixType sparseLu;
    sparseLu = spLtwisted.block(b_rows, b_rows, b_cols, b_cols);

    Eigen::VectorXd xs(b_rows);
    xs.setConstant(1);
    xs.segment(0, idxBackground.size()).setConstant(0);

    SparseMatrixType sparseLHard;
    sparseLHard = sparseBt * sparseB + sparseLu * sparseLu;
    Eigen::VectorXd bHard = -(sparseLu * sparseBt + sparseBt * sparseLs) * xs;

    Eigen::CholmodDecomposition<SparseMatrixType, Eigen::Lower> sparseSolver;
    sparseSolver.setMode(Eigen::CholmodSupernodalLLt);
    sparseSolver.compute(sparseLHard);

    m_solutionVector = sparseSolver.solve(bHard);
}
