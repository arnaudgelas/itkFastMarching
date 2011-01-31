/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009-11

 Copyright (c) 2009-11, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#ifndef READPATHFILE_H
#define READPATHFILE_H

/////////////////////////////////////////////////////////////
// Reads a *.path file and adds the path info to the given filter
template <class PathFilterType>
int ReadPathFile( const char * PathFilename, typename PathFilterType::Pointer pathFilter )
{
    // Path file example:
    // Path: [272.00, 128.00] [490.00, 148.00]
    // Path: [272.00, 128.00] [381.00, 001.00]
    // Path: [272.00, 128.00] [002.00, 130.00]
    // Path: [272.00, 128.00] [274.00, 268.00]

    // NOTE: No checking is done on the path file: the user must ensure it is valid!!!
    std::string filename = PathFilename;
    if ( !itksys::SystemTools::FileIsFullPath(PathFilename) )
      {
      std::string currentDirectory = itksys::SystemTools::GetCurrentWorkingDirectory();
      filename = currentDirectory + "/" + filename;
      }

    std::ifstream file(filename.c_str(), std::ios::in);
    if ( !file )
      {
      std::cerr << "Unable to open file ";
      std::cerr << PathFilename;
      std::cerr << " for reading.";
      std::cerr << std::endl;
      return EXIT_FAILURE;
      }

    std::string line;
    bool has_newline = false;
    while ( itksys::SystemTools::GetLineFromStream(file, line, &has_newline) )
      {
      if (has_newline)
        {
        typename PathFilterType::PathInfo info;
        itksys::SystemTools::ReplaceString( line, "Path: ", "" );
        itksys::SystemTools::ReplaceString( line, " ", "" );
        itksys::SystemTools::ReplaceString( line, "[", "" );
        std::vector<itksys::String> parts;
        parts = itksys::SystemTools::SplitString( line.c_str(), ']' );
        unsigned int numNonNullParts = 0;
        for (unsigned int i=0; i<parts.size(); i++)
          {
          if ( parts[i].length() != 0 )
            {
            numNonNullParts++;
            }
          }

        for (unsigned int i=0; i<numNonNullParts; i++)
          {
          if ( parts[i].length() != 0 )
            {
            typename PathFilterType::PointType point;
            std::vector<itksys::String> partsPoint;
            partsPoint = itksys::SystemTools::SplitString( parts[i].c_str(), ',' );
            for (unsigned int j=0; j<partsPoint.size(); j++)
              {
              point[j] = atof( partsPoint[j].c_str() );
              }
            if ( i==0 )
              {
              info.SetStartPoint( point );
              }
            else
              {
              if ( i == numNonNullParts - 1 )
                {
                info.SetEndPoint( point );
                }
              else
                {
                info.AddWayPoint( point );
                }
              }
            }
          }
        pathFilter->AddPathInfo( info );
        }
      }

    return EXIT_SUCCESS;
}

#endif // READPATHFILE_H
