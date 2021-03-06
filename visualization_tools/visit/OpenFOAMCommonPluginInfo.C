/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#include <OpenFOAMPluginInfo.h>
#include <avtOpenFOAMFileFormat.h>
#include <avtMTMDFileFormatInterface.h>
#include <avtGenericDatabase.h>

// ****************************************************************************
//  Method:  OpenFOAMCommonPluginInfo::GetDatabaseType
//
//  Purpose:
//    Returns the type of a OpenFOAM database.
//
//  Programmer:  root -- generated by xml2info
//  Creation:    Mon Jul 24 10:00:07 PDT 2006
//
// ****************************************************************************
DatabaseType
OpenFOAMCommonPluginInfo::GetDatabaseType()
{
    return DB_TYPE_MTMD;
}

// ****************************************************************************
//  Method:  OpenFOAMCommonPluginInfo::GetDefaultExtensions
//
//  Purpose:
//    Returns the default extensions for a OpenFOAM database.
//
//  Programmer:  root -- generated by xml2info
//  Creation:    Mon Jul 24 10:00:07 PDT 2006
//
// ****************************************************************************
std::vector<std::string>
OpenFOAMCommonPluginInfo::GetDefaultExtensions()
{
    std::vector<std::string> defaultExtensions;

    return defaultExtensions;
}

// ****************************************************************************
//  Method: OpenFOAMCommonPluginInfo::SetupOpenFOAMDatabase
//
//  Purpose:
//      Sets up a OpenFOAM database.
//
//  Arguments:
//      list    A list of file names.
//      nList   The number of timesteps in list.
//      nBlocks The number of blocks in the list.
//
//  Returns:    A OpenFOAM database from list.
//
//  Programmer: root -- generated by xml2info
//  Creation:   Mon Jul 24 10:00:07 PDT 2006
//
// ****************************************************************************
avtDatabase *
OpenFOAMCommonPluginInfo::SetupDatabase(const char *const *list,
                                   int nList, int nBlock)
{
    return new avtGenericDatabase(
               new avtMTMDFileFormatInterface(
                   new avtOpenFOAMFileFormat(list[0])));
}

// ****************************************************************************
//  Method: OpenFoamCommonPluginInfo::GetFilenames
//
//  Purpose:
//      Allows the user to open files that have no extensions
//
//  Returns:    A vector of filenames
//
//  Programmer: root -- generated by xml2info
//  Creation:   Mon Jul 24 10:00:07 PDT 2006
//
// ****************************************************************************
std::vector<std::string>
OpenFOAMCommonPluginInfo::GetFilenames()
{
   std::vector<std::string> filenames;
   filenames.push_back("controlDict");

   return filenames;
}
