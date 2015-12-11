/***************************************************************************
 *  Project:    netpop
 *  File:       utils.c++
 *  Language:   C++
 *
 *  netpop is free software: you can redistribute it and/or modify it under the
 *  terms of the GNU General Public License as published by the Free Software
 *  Foundation, either version 3 of the License, or (at your option) any later
 *  version.
 *
 *  netpop is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public License along with
 *  NeutralClusters.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright   Mark Padgham December 2015
 *  Author:     Mark Padgham
 *  E-Mail:     mark.padgham@email.com
 *
 *  Description:    Simulates populations reproducing according to cubic
 *                  logistic dynamic and moving throughout weighted and directed
 *                  networks of three and four nodes.
 *
 *  Project Structure:  
 *      Routines are divided between the two main programs:
 *      1. net3:    Triangular network of 3 nodes
 *      2. net4:    Four kinds of networks of 4 nodes (see net.c++ for
 *                  diagrams).
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options
 ***************************************************************************/


#include "utils.h"

void progLine (double progress, double tseconds)
{
    // progress bar with progress in [0,1] and prefix integer [nc] and elapsed
    // time [tseconds]
    struct winsize w;
    ioctl (0, TIOCGWINSZ, &w); 
    int ncols = w.ws_col; // Number of columns in console
    ncols -= 10;

    int barlen = ncols - 20;
    int proglen = floor (barlen * progress);
    int gaplen = barlen - proglen;

    std::string time_taken = get_time (tseconds);
    std::cout << "[" << time_taken << "] |";
    for (int i=0; i<proglen; i++) 
        std::cout << "-";
    for (int i=0; i<gaplen; i++) 
        std::cout << " ";
    std::cout << "| " << (int) floor (progress * 100.0) << "%\r";
    std::cout.flush ();
}

std::string get_time (double tseconds)
{
    std::stringstream sstr;
    sstr.str ("");
    
    int hh = floor (tseconds / 3600.0);
    if (hh == 0) 
        sstr << "00:";
    else if (hh < 10) 
        sstr << "0" << hh << ":";
    else 
        sstr << hh << ":";
    double trem = tseconds - (double) hh * 3600.0;
    int mm = floor (trem / 60.0);
    if (mm == 0) 
        sstr << "00:";
    else if (mm < 10) 
        sstr << "0" << mm << ":";
    else 
        sstr << mm << ":";
    int ss = round (trem - (double) mm * 60.0);
    if (ss == 0) 
        sstr << "00";
    else if (ss < 10) 
        sstr << "0" << ss;
    else 
        sstr << ss;

    return (sstr.str ());
} // end function timeout
