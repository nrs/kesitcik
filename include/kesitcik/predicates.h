/* GTS - Library for the manipulation of triangulated surfaces
 * Copyright (C) 1999 St�phane Popinet
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
/* Header file for robust predicates by Jonathan Richard Shewchuk */

#ifndef __PREDICATES_H__
#define __PREDICATES_H__

#define REAL double                      /* float or double */


namespace predicates 
{
REAL exactinit();

REAL orient2d (REAL * pa,REAL * pb,REAL * pc);
REAL orient2dexact (REAL * pa,REAL * pb,REAL * pc);
//REAL orient2d            (REAL * pa,REAL * pb,REAL * pc);
REAL orient3d (REAL * pa,REAL * pb,REAL * pc,
			    REAL * pd);
REAL incircle (REAL * pa,REAL * pb,REAL * pc,
			    REAL * pd);
REAL insphere (REAL * pa,REAL * pb,REAL * pc,REAL * pd,REAL * pe);
//REAL orient2d(REAL pa[2],REAL pb[2],REAL pc[2]);
}
#endif /* __PREDICATES_H__ */
