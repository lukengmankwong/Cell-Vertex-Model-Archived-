#ifndef LIBRARIES_H
#define LIBRARIES_H

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel                   K;
typedef CGAL::Vector_2<K>                                                   Vec;
typedef CGAL::Point_2<K>												  Point;

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;

#endif // LIBRARIES_H
