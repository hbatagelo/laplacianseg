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

#ifndef COMMANDLINEPARSER_H
#define COMMANDLINEPARSER_H

#include <QColor>
#include <QCommandLineParser>

class CommandLineParser
{
public:
    enum Status
    {
        CommandLineOk,
        CommandLineError,
        CommandLineVersionRequested,
        CommandLineHelpRequested
    };

private:
    QCommandLineParser m_parser;

    QString m_inputImagePath;
    QString m_inputSeedsPath;
    QString m_outputPath;

    bool m_useHardLabels;
    bool m_useSuperpixels;
    bool m_outputAsBinary;
    bool m_quietMode;
    int m_spSize;
    double m_compactness;

    QColor m_foregroundSeedColor = QColor(Qt::red);
    QColor m_backgroundSeedColor = QColor(Qt::blue);

    Status m_status = CommandLineOk;
    QString m_statusMessage;

public:
    CommandLineParser();

    QString inputImagePath() const { return m_inputImagePath; }
    QString inputSeedsPath() const { return m_inputSeedsPath; }
    QString outputPath() const { return m_outputPath; }
    bool useHardLabels() const { return m_useHardLabels; }
    bool useSuperpixels() const { return m_useSuperpixels; }
    bool outputAsBinary() const { return m_outputAsBinary; }
    bool quietMode() const { return m_quietMode; }
    QColor foregroundSeedColor() const { return  m_foregroundSeedColor; }
    QColor backgroundSeedColor() const { return  m_backgroundSeedColor; }
    int superpixelSize() const { return m_spSize; }
    double compactness() const { return m_compactness; }

    Status status() const { return m_status; }
    QString statusMessage() const { return m_statusMessage; }

    QString helpText() const { return m_parser.helpText(); }
    [[noreturn]] void showHelp(int exitCode = 0) { m_parser.showHelp(exitCode); }
};

#endif // COMMANDLINEPARSER_H
