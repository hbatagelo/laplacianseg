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

#include "commandlineparser.h"

CommandLineParser::CommandLineParser()
{
    m_parser.setApplicationDescription(
        "Seeded image segmentation using Laplacian Coordinates.\n"
        "(c) 2020 Harlen Batagelo, João Paulo Gois");

    m_parser.setSingleDashWordOptionMode(QCommandLineParser::ParseAsLongOptions);

    m_parser.addPositionalArgument("input", "Input image file.");
    m_parser.addPositionalArgument("seeds", "Input seeds image file.");
    m_parser.addPositionalArgument("output", "Output image file.");

    QCommandLineOption fgColorOption(
        "fg", "Sets color of foreground seeds (default is '#ff0000').", "color", "#ff0000");
    m_parser.addOption(fgColorOption);

    QCommandLineOption bgColorOption(
        "bg", "Sets color of background seeds (default is '#0000ff').", "color", "#0000ff");
    m_parser.addOption(bgColorOption);

    QCommandLineOption hardOption(
        "hard", "Uses seeds as hard labeling constraints.");
    m_parser.addOption(hardOption);

    QCommandLineOption superpixelOption(
        "superpixel", "Uses SLIC superpixels.");
    m_parser.addOption(superpixelOption);

    QCommandLineOption spSizeOption(
        "size", "Sets superpixel size (default is 100).", "value", "100");
    m_parser.addOption(spSizeOption);

    QCommandLineOption spCompactnessOption(
        "compactness", "Sets superpixel compactness (default is 10.0).", "value", "10.0");
    m_parser.addOption(spCompactnessOption);

    QCommandLineOption binaryOption(
        {"b", "binary"}, "Writes output as a binary image.");
    m_parser.addOption(binaryOption);

    QCommandLineOption quietOption(
        {"q", "quiet"}, "Runs in silent mode.");
    m_parser.addOption(quietOption);

    QCommandLineOption versionOption = m_parser.addVersionOption();

    if (!m_parser.parse(QCoreApplication::arguments()))
    {
        m_statusMessage = m_parser.errorText();
        m_status = CommandLineError;
        return;
    }

    if (m_parser.isSet(versionOption))
    {
        m_status = CommandLineVersionRequested;
        return;
    }

    QStringList positionalArguments = m_parser.positionalArguments();
    if (positionalArguments.isEmpty())
    {
        m_statusMessage = "Argument 'input' missing.";
        m_status = CommandLineError;
        return;
    }
    else if (positionalArguments.size() == 1)
    {
        m_statusMessage = "Argument 'seeds' missing.";
        m_status = CommandLineError;
        return;
    }
    else if (positionalArguments.size() == 2)
    {
        m_statusMessage = "Argument 'output' missing.";
        m_status = CommandLineError;
        return;
    }

    QStringList args = m_parser.positionalArguments();
    m_inputImagePath = args.at(0);
    m_inputSeedsPath = args.at(1);
    m_outputPath = args.at(2);
    m_useHardLabels = m_parser.isSet(hardOption);
    m_useSuperpixels = m_parser.isSet(superpixelOption);
    m_outputAsBinary = m_parser.isSet(binaryOption);
    m_quietMode = m_parser.isSet(quietOption);
    m_spSize = m_parser.value(spSizeOption).toInt();
    if (m_spSize < 1)
    {
        m_statusMessage = "Invalid superpixel size.";
        m_status = CommandLineError;
        return;
    }
    m_compactness = m_parser.value(spCompactnessOption).toDouble();
    if (m_compactness < 1.0)
    {
        m_statusMessage = "Invalid compactness.";
        m_status = CommandLineError;
        return;
    }

    QString fgColorName = m_parser.value(fgColorOption);
    if (QColor::isValidColor(fgColorName))
    {
        m_foregroundSeedColor = QColor(fgColorName);
    }
    else
    {
        m_statusMessage = "Invalid foreground seed color.";
        m_status = CommandLineError;
        return;
    }

    QString bgColorName = m_parser.value(bgColorOption);
    if (QColor::isValidColor(bgColorName))
    {
        m_backgroundSeedColor = QColor(bgColorName);
    }
    else
    {
        m_statusMessage = "Invalid background seed color.";
        m_status = CommandLineError;
        return;
    }
}
