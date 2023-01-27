#include<surface_integration.h>
#include<sstream>

namespace surfaceIntegration{

    errorOut decomposeSphere( const floatType &radius, unsigned int &elementCount,
                              floatVector &points, std::vector<unsigned int> &connectivity ){
        /*!
         * Decompose a sphere into quadratic elements for use in numeric integration
         * 
         * The decomposition takes place by projecting a unit cube onto a sphere of the given radius
         * the element count is the number of elements on the side of the cube's face. The nodes are
         * defined following the convention detailed in "The Finite Element Method," T. J. R. Hughes
         * Chapter 3, Section 6: Higher Order Elements; Lagrange Polynomials
         * 
         * [4] - [7] - [3]
         *  |     |     |
         * [8] - [9] - [6]
         *  |     |     |
         * [1] - [5] - [2]
         * 
         * \param &radius: The radius of the sphere
         * \param &elementCount: The number of elements on each edge of the projected unit cube
         * \param &points: The points on the surface of the sphere (p1x, p1y, p1z, p2x, p2y, p2z, ...)
         * \param &connectivity: The connectivity array for the elements (e1_1, e1_2, e1_3, e1_4, e1_5, e1_6, e1_7, e1_8, e1_9, e2_1, ...
         */

        // Set the number of points on the unit cube face edges
        unsigned int n_points_edge = 2 * elementCount + 1;

        return NULL;

    }

    errorOut buildSurfacePoints( const floatType &x0, const floatType &y0, const floatType &z0,
                                 const floatType &dx, const floatType &dy,
                                 const unsigned int n_points_x, const unsigned int n_points_y,
                                 floatVector &points ){
        /*!
         * Build a collection of points which represent a planar surface
         * 
         * \param &x0: The starting point in the x direction
         * \param &y0: The starting point in the y direction
         * \param &z0: The z position (constant for all points)
         * \param &dx: The spacing between points in the x direction
         * \param &dy: The spacing between points in the y direction
         * \param &n_points_x: The number of points in the x direction
         * \param &n_points_y: The number of points in the y direction
         * \param &points: The resulting surface points
         */

        floatType x = x0;
        floatType y = y0;
        floatType z = z0;
        
        points = floatVector( 3 * n_points_x * n_points_y, 0 );
        
        unsigned int index = 0;
        
        for ( unsigned int i = 0; i < n_points_y; i++ ){

            for ( unsigned int j = 0; j < n_points_x; j++ ){

                points[index + 0] = x;
                points[index + 1] = y;
                points[index + 2] = z;
                
                index += 3;
                
                x += dx;

            }

            x = x0;

            y += dy;

        }

        return NULL;

    }

    errorOut rotatePoints( const floatVector &points,
                           const floatType &thetaX, const floatType &thetaY, const floatType &thetaZ,
                           floatVector &rotatedPoints ){
        /*!
         * Rotate the vector of points in row-major vector form to a new configuration
         * 
         * \f$x_i' = R_{ij}^z R_{jk}^y R_{kl}^x x_l\f$
         * 
         * \param &points: The points in the initial configuration
         * \param &thetaX: The rotation about the X axis
         * \param &thetaY: The rotation about the Y axis
         * \param &thetaZ: The rotation about the Z axis
         * \param &rotatedPoints: The rotated points
         */

        unsigned int n_points = points.size( ) / 3;

        if ( n_points * 3 != points.size( ) ){

            return new errorNode( __func__, "The size of points is not a multiple of 3:\n  points.size( ) = " + std::to_string( points.size( ) ) );

        }

        floatType cx = std::cos( thetaX );
        floatType sx = std::sin( thetaX );

        floatType cy = std::cos( thetaY );
        floatType sy = std::sin( thetaY );

        floatType cz = std::cos( thetaZ );
        floatType sz = std::sin( thetaZ );
        
        floatMatrix Rx = { 
                             { 1,  0,   0 },
                             { 0, cx, -sx },
                             { 0, sx,  cx },
                         };
         
        floatMatrix Ry = {
                             { cy, 0, sy },
                             {  0, 1,  0 },
                             {-sy, 0, cy },
                         };
         
        floatMatrix Rz = {
                             { cz, -sz, 0 },
                             { sz,  cz, 0 },
                             {  0,   0, 1 }
                         };
         
        floatMatrix R = vectorTools::dot( Rz, vectorTools::dot( Ry, Rx ) );

        rotatedPoints = floatVector( points.size( ), 0 );

        floatVector temp_vector;

        for ( unsigned int index = 0; index < n_points; index++ ){

            temp_vector = vectorTools::dot( R, floatVector( points.begin( ) + 3 * index, points.begin( ) + 3 * ( index + 1 ) ) );

            rotatedPoints[ 3 * index + 0 ] = temp_vector[ 0 ];
            rotatedPoints[ 3 * index + 1 ] = temp_vector[ 1 ];
            rotatedPoints[ 3 * index + 2 ] = temp_vector[ 2 ];

        }

        return NULL;

    }

}
