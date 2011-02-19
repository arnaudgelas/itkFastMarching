/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

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
