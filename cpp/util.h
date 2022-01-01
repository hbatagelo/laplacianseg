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

#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include <concepts>
#include <numeric>
#include <opencv2/opencv.hpp>
#include <ranges>
#include <span>

namespace lc {

// In-place RGB to Lab
template <std::floating_point T> void RGBToLab(std::span<T, 3> RGB) {
  // Inverse sRGB companding
  std::array<T, 3> linRGB{};
  std::ranges::transform(RGB, linRGB.begin(), [](T t) {
    return t > 0.04045F ? std::pow((t + 0.055F) / 1.055F, 2.4F) : t / 12.92F;
  });

  // Convert RGB to X (idx=0), Y (idx=1) or Z (idx=2)
  auto RGBToXYZ{[](auto const &rgb, std::size_t const idx) {
    // Illuminant = D65 (2°, CIE 1931)
    static std::array const refWhite{0.95047F, 1.0F, 1.08883F};
    // RGB->XYZ conversion matrix (sRGB - D65)
    // clang-format off
    static std::array const RGBToXYZ{
      0.4124564F, 0.3575761F, 0.1804375F,
      0.2126729F, 0.7151522F, 0.0721750F,
      0.0193339F, 0.1191920F, 0.9503041F};
    // clang-format on
    auto const &s{std::span{RGBToXYZ}.subspan(idx * 3, 3)};
    return std::inner_product(s.begin(), s.end(), rgb.begin(), T{}) /
           refWhite.at(idx);
  }};
  std::array XYZ{RGBToXYZ(linRGB, 0), RGBToXYZ(linRGB, 1), RGBToXYZ(linRGB, 2)};

  // XYZ->Lab
  std::ranges::transform(XYZ, XYZ.begin(), [](T t) {
    return t > 0.008856F ? std::pow(t, 1 / 3.0F) : 7.787F * t + (16 / 116.0F);
  });
  auto const &[fX, fY, fZ]{XYZ};
  RGB[0] = 116.0F * fY - 16.0F;
  RGB[1] = 500.0F * (fX - fY);
  RGB[2] = 200.0F * (fY - fZ);
}

template <std::floating_point T>
constexpr bool almostEquals(T const &a, T const &b) {
  return std::abs(a - b) <= (std::max(std::abs(a), std::abs(b)) *
                             std::numeric_limits<T>::epsilon());
}

template <typename T>
concept xy_coordinates = requires(T a) {
  std::is_arithmetic<decltype(a.x())>();
  std::is_arithmetic<decltype(a.y())>();
};

template <xy_coordinates T>
constexpr bool inBounds(T const &p, T const &lo, T const &hi) {
  return p.x() >= lo.x() && p.x() < hi.x() && p.y() >= lo.y() && p.y() < hi.y();
}

// Apply image thresholding to compute the binary image. If otsuVariant is true,
// the threshold is given by Isol > (1-otsu) instead of Isol > otsu
inline cv::Mat ApplyThresholding(cv::Mat const &image, bool otsu = true,
                                 bool otsuVariant = true) {
  auto const binary{cv::THRESH_BINARY};

  cv::Mat output(image.size(), image.type());

  if (otsu) {
    auto const binaryOtsu{binary | cv::THRESH_OTSU};
    if (otsuVariant) {
      auto otsuThreshold{cv::threshold(image, output, 0, 255, binaryOtsu)};
      if (almostEquals(otsuThreshold, 0.0))
        otsuThreshold = 255;
      cv::threshold(image, output, 255 - otsuThreshold, 255, binary);
    } else {
      cv::threshold(image, output, 0, 255, binaryOtsu);
    }
  } else {
    cv::threshold(image, output, 127, 255, binary);
  }
  return output;
}

} // namespace lc

#endif // UTIL_H
