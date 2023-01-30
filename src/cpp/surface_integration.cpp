#include<surface_integration.h>
#include<sstream>

namespace surfaceIntegration{

    errorOut decomposeSphere( const floatType &radius, const unsigned int &elementCount,
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

        errorOut error;

        // Form the base structure of the sphere
        floatVector cubePoints;
        error = formBaseCubePoints( elementCount, cubePoints );
        if ( error ){
            errorOut result = new errorNode( __func__, "Error in the formation of the base cube" );
            result->addNext( error );
            return error;
        }

        // Re-scale the points to be a sphere
        points = floatVector( cubePoints.size( ), 0 );
        unsigned int n_points = cubePoints.size( ) / 3;

        if ( n_points * 3 != cubePoints.size( ) ){
            return new errorNode( __func__, "The length of cube points isn't consistent with a collection of 3D points\n  cubePoints.size( ): " + std::to_string( cubePoints.size( ) ) );
        }

        for ( unsigned int i = 0; i < n_points; i++ ){

            floatType length = vectorTools::l2norm( floatVector( cubePoints.begin( ) + 3 * i, cubePoints.begin( ) + 3 * ( i + 1 ) ) );

            points[ 3 * i + 0 ] = radius * cubePoints[ 3 * i + 0 ] / length;
            points[ 3 * i + 1 ] = radius * cubePoints[ 3 * i + 1 ] / length;
            points[ 3 * i + 2 ] = radius * cubePoints[ 3 * i + 2 ] / length;

        }

        // From the connectivity
        error = formCubeConnectivity( elementCount, connectivity );
        if ( error ){
            errorOut result = new errorNode( __func__, "Error in the formation of the base cube connectivity" );
            result->addNext( error );
            return error;
        }

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

    errorOut formSurfaceConnectivity( const std::vector< unsigned int > &surfaceIDs,
                                      const unsigned int &n_elements_x, const unsigned int &n_elements_y,
                                      unsigned int &index, std::vector< unsigned int > &connectivity ){
        /*!
         * Form the connectivity for a planar surface of quadratic elements
         * 
         * \param &surfaceIDs: The id numbers of the surface nodes
         * \param &n_elements_x: The number of elements in the local x direction
         * \param &n_elements_y: The number of elements in the local y direction
         * \param &index: The most recent element number
         * \param &connectivity: The connectivity array
         */

        if ( connectivity.size( ) < ( 9 * n_elements_x * n_elements_y + 9 * index ) ){

            return new errorNode( __func__, "The connectivity array is too small for the expected number of elements" );

        }

        unsigned int delta_p = 3; // The number of nodes on each side of the element
                                
        for ( unsigned int j = 0; j < n_elements_y; j++ ){

            for ( unsigned int i = 0; i < n_elements_x; i++ ){
                
                connectivity[ 9 * index + 0 ] = surfaceIDs[ 2 * ( 2 * n_elements_x + 1 ) * j + ( delta_p - 1 ) * i ];
                connectivity[ 9 * index + 1 ] = surfaceIDs[ 2 * ( 2 * n_elements_x + 1 ) * j + ( delta_p - 1 ) * i + 2 ];
                connectivity[ 9 * index + 2 ] = surfaceIDs[ 2 * ( 2 * n_elements_x + 1 ) * j + ( delta_p - 1 ) * i + 2 * ( 2 * n_elements_x + 1 ) + 2 ];
                connectivity[ 9 * index + 3 ] = surfaceIDs[ 2 * ( 2 * n_elements_x + 1 ) * j + ( delta_p - 1 ) * i + 2 * ( 2 * n_elements_x + 1 ) ];
                
                connectivity[ 9 * index + 4 ] = surfaceIDs[ 2 * ( 2 * n_elements_x + 1 ) * j + ( delta_p - 1 ) * i + 1 ];
                connectivity[ 9 * index + 5 ] = surfaceIDs[ 2 * ( 2 * n_elements_x + 1 ) * j + ( delta_p - 1 ) * i + ( 2 * n_elements_x + 1 ) + 2 ];
                connectivity[ 9 * index + 6 ] = surfaceIDs[ 2 * ( 2 * n_elements_x + 1 ) * j + ( delta_p - 1 ) * i + 2 * ( 2 * n_elements_x + 1 ) + 1 ];
                connectivity[ 9 * index + 7 ] = surfaceIDs[ 2 * ( 2 * n_elements_x + 1 ) * j + ( delta_p - 1 ) * i + ( 2 * n_elements_x + 1 ) ];
                
                connectivity[ 9 * index + 8 ] = surfaceIDs[ 2 * ( 2 * n_elements_x + 1 ) * j + ( delta_p - 1 ) * i + ( 2 * n_elements_x + 1 ) + 1 ];
                
                index += 1;

            }

        }

        return NULL;

    }

    errorOut formBaseCubePoints( const unsigned int &elementCount, floatVector &points ){
        /*!
         * Form the base cube points which will be used for integration using quadratic elements
         * 
         * \param &elementCount: The number of elements on each edge of the cube
         * \param &points: The resulting points
         */

        errorOut error = NULL;
        floatType pi = 3.141592653589793;
        unsigned int n_points_edge = 2 * elementCount + 1;

        floatType x = -1;
        floatType y = -1;
        floatType z =  1;

        floatType dx = 1. / elementCount;
        floatType dy = 1. / elementCount;

        // Build the top surface
        floatVector top_points;
        error = buildSurfacePoints( x, y, z, dx, dy, n_points_edge, n_points_edge - 1, top_points );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when building the top surface" );
            result->addNext( error );
            return result;

        }

        // Build the back surface
        floatVector temp;
        floatVector back_points;
        error = buildSurfacePoints( x, y, z, dx, dy, n_points_edge, n_points_edge - 1, temp );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when building the back surface" );
            result->addNext( error );
            return result;

        }

        error = rotatePoints( temp, -0.5 * pi, 0, 0, back_points );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when rotating the back surface" );
            result->addNext( error );
            return result;

        }

        // Build the bottom surface
        floatVector bottom_points;
        error = buildSurfacePoints( x, y, z, dx, dy, n_points_edge, n_points_edge - 1, temp );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when building the bottom surface" );
            result->addNext( error );
            return result;

        }

        error = rotatePoints( temp, -pi, 0, 0, bottom_points );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when rotating the bottom surface" );
            result->addNext( error );
            return result;

        }

        // Build the front surface
        floatVector front_points;
        error = buildSurfacePoints( x, y, z, dx, dy, n_points_edge, n_points_edge - 1, temp );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when building the front surface" );
            result->addNext( error );
            return result;

        }

        error = rotatePoints( temp, -1.5 * pi, 0, 0, front_points );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when rotating the front surface" );
            result->addNext( error );
            return result;

        }

        // Build the right surface
        floatVector right_points;
        error = buildSurfacePoints( x + dx, y + dy, z, dx, dy, n_points_edge - 2, n_points_edge - 2, temp );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when building the right surface" );
            result->addNext( error );
            return result;

        }

        error = rotatePoints( temp, 0, 0.5 * pi, 0, right_points );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when rotating the right surface" );
            result->addNext( error );
            return result;

        }

        // Build the left surface
        floatVector left_points;
        error = buildSurfacePoints( x + dx, y + dy, z, dx, dy, n_points_edge - 2, n_points_edge - 2, temp );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when building the left surface" );
            result->addNext( error );
            return result;

        }

        error = rotatePoints( temp, 0, -0.5 * pi, 0, left_points );
        if ( error ){

            errorOut result = new errorNode( __func__, "Error when rotating the left surface" );
            result->addNext( error );
            return result;

        }

        points = vectorTools::appendVectors( { top_points, back_points, bottom_points, front_points, right_points, left_points } );

        return NULL;

    }

    errorOut formCubeConnectivity( const unsigned int &elementCount, std::vector< unsigned int > &connectivity ){
        /*!
         * Form the connectivity vector for the base cube defined through formBaseCubePoints
         * 
         * \param &elementCount: The number of elements along each edge of the base cube
         * \param &connectivity: The connectivity vector for the quadratic elements
         */

        connectivity = std::vector< unsigned int >( 9 * 6 * elementCount * elementCount );
        errorOut error;

        unsigned int n_points_edge = 2 * elementCount + 1;

        unsigned int index = 0;

        // Define the top surface connectivity
        std::vector< unsigned int > surfaceIDs( n_points_edge * n_points_edge );
        std::iota( surfaceIDs.begin( ), surfaceIDs.end( ), 0 );
        error = formSurfaceConnectivity( surfaceIDs, elementCount, elementCount, index, connectivity );
        if ( error ){
            errorOut result = new errorNode( __func__, "Error when building the connectivity of the top surface" );
            result->addNext( error );
            return result;
        }

        // Define the back surface connectivity
        std::iota( surfaceIDs.begin( ), surfaceIDs.end( ), n_points_edge * ( n_points_edge - 1 ) );
        error = formSurfaceConnectivity( surfaceIDs, elementCount, elementCount, index, connectivity );
        if ( error ){
            errorOut result = new errorNode( __func__, "Error when building the connectivity of the top surface" );
            result->addNext( error );
            return result;
        }

        // Define the bottom surface connectivity
        std::iota( surfaceIDs.begin( ), surfaceIDs.end( ), 2 * n_points_edge * ( n_points_edge - 1 ) );
        error = formSurfaceConnectivity( surfaceIDs, elementCount, elementCount, index, connectivity );
        if ( error ){
            errorOut result = new errorNode( __func__, "Error when building the connectivity of the top surface" );
            result->addNext( error );
            return result;
        }

        // Define the front surface connectivity
        std::iota( surfaceIDs.begin( ), surfaceIDs.end( ) - n_points_edge, 3 * n_points_edge * ( n_points_edge - 1 ) );
        std::iota( surfaceIDs.end( ) - n_points_edge, surfaceIDs.end( ), 0 );
        error = formSurfaceConnectivity( surfaceIDs, elementCount, elementCount, index, connectivity );
        if ( error ){
            errorOut result = new errorNode( __func__, "Error when building the connectivity of the top surface" );
            result->addNext( error );
            return result;
        }

        // Define the right surface connectivity
        unsigned int offset = 4 * n_points_edge * ( n_points_edge - 1 ) - 1;
        std::vector< unsigned int > bottom_edge( n_points_edge );
        bottom_edge[ 0 ] = n_points_edge - 1;
        for ( unsigned int i = 1; i < bottom_edge.size( ); i++ ){
            bottom_edge[ i ] = offset - n_points_edge * ( i - 1 );
        }

        offset = n_points_edge * ( n_points_edge - 1 );
        std::vector< unsigned int > top_edge( n_points_edge );
        for ( unsigned int i = 0; i < top_edge.size( ); i++ ){
            top_edge[ i ] = offset + n_points_edge * i + ( n_points_edge - 1 );
        }

        std::vector< unsigned int > left_edge( n_points_edge - 2 );
        for ( unsigned int i = 0; i < left_edge.size( ); i++ ){
            left_edge[ i ] = n_points_edge * ( i + 1 ) + n_points_edge - 1;
        }

        std::vector< unsigned int > right_edge( n_points_edge - 2 );
        offset = 3 * n_points_edge * ( n_points_edge - 1 ) - 1;
        for ( unsigned int i = 0; i < left_edge.size( ); i++ ){
            right_edge[ i ] = offset - n_points_edge * i;
        }

        offset = 4 * n_points_edge * ( n_points_edge - 1 );
        std::vector< unsigned int > center( ( n_points_edge - 2 ) * ( n_points_edge - 2 ) );
        std::iota( center.begin( ), center.end( ), offset );

        for ( unsigned int i = 0; i < n_points_edge; i++ ){
            surfaceIDs[ i ] = bottom_edge[ i ];
            surfaceIDs[ n_points_edge * ( n_points_edge - 1 ) + i ] = top_edge[ i ];
        }

        for ( unsigned int i = 0; i < ( n_points_edge - 2 ); i++ ){

            surfaceIDs[ n_points_edge * i + bottom_edge.size( ) ] = left_edge[ i ];

            for ( unsigned int j = 0; j < ( n_points_edge - 2 ); j++ ){

                surfaceIDs[ n_points_edge * i + bottom_edge.size( ) + j + 1 ] = center[ ( n_points_edge - 2 ) * i + j ];

            }

            surfaceIDs[ n_points_edge * i + bottom_edge.size( ) + n_points_edge - 1 ] = right_edge[ i ];

        }
        error = formSurfaceConnectivity( surfaceIDs, elementCount, elementCount, index, connectivity );
        if ( error ){
            errorOut result = new errorNode( __func__, "Error when building the connectivity of the right surface" );
            result->addNext( error );
            return result;
        }

        // Define the left surface connectivity
        offset = 3 * n_points_edge * ( n_points_edge - 1 );
        for ( unsigned int i = 0; i < ( bottom_edge.size( ) - 1 ); i++ ){
            bottom_edge[ i ] = offset + n_points_edge * i;
        }
        bottom_edge[ bottom_edge.size( ) - 1 ] = 0;

        offset = 2 * n_points_edge * ( n_points_edge - 1 );
        for ( unsigned int i = 0; i < top_edge.size( ); i++ ){
            top_edge[ i ] = offset  - n_points_edge * i;
        }
        
        for ( unsigned int i = 0; i < ( n_points_edge - 2 ); i++ ){
            right_edge[ i ] = n_points_edge * ( i + 1 );
        }

        offset = 3 * n_points_edge * ( n_points_edge - 1 );
        for ( unsigned int i = 0; i < ( n_points_edge - 2 ); i++ ){
            left_edge[ i ] = offset - n_points_edge * ( i + 1 );
        }

        offset = 4 * n_points_edge * ( n_points_edge - 1 ) + ( n_points_edge - 2 ) * ( n_points_edge - 2 );
        std::iota( center.begin( ), center.end( ), offset );

        for ( unsigned int i = 0; i < n_points_edge; i++ ){
            surfaceIDs[ i ] = bottom_edge[ i ];
            surfaceIDs[ n_points_edge * ( n_points_edge - 1 ) + i ] = top_edge[ i ];
        }

        for ( unsigned int i = 0; i < ( n_points_edge - 2 ); i++ ){

            surfaceIDs[ n_points_edge * i + bottom_edge.size( ) ] = left_edge[ i ];

            for ( unsigned int j = 0; j < ( n_points_edge - 2 ); j++ ){

                surfaceIDs[ n_points_edge * i + bottom_edge.size( ) + j + 1 ] = center[ ( n_points_edge - 2 ) * i + j ];

            }

            surfaceIDs[ n_points_edge * i + bottom_edge.size( ) + n_points_edge - 1 ] = right_edge[ i ];

        }
        error = formSurfaceConnectivity( surfaceIDs, elementCount, elementCount, index, connectivity );
        if ( error ){
            errorOut result = new errorNode( __func__, "Error when building the connectivity of the left surface" );
            result->addNext( error );
            return result;
        }

        return NULL;

    }

    errorOut evaluateQuadraticShapeFunctions( const floatType &xi, const floatType &eta, floatVector &shapeFunctions ){
        /*!
         * Evaluate the shape functions of a quadratic element
         * 
         * \param &xi: Local xi coordinate
         * \param &eta: Local eta coordinate
         * \param &shapeFunctions: The resulting shape functions
         */

        floatMatrix localPoints = {
                                      { -1, -1 },
                                      {  1, -1 },
                                      {  1,  1 },
                                      { -1,  1 },
                                      {  0, -1 },
                                      {  1,  0 },
                                      {  0,  1 },
                                      { -1,  0 },
                                      {  0,  0 }
                                  };

        shapeFunctions = floatVector( 9, 0 );

        // Corner nodes
        for ( unsigned int i = 0; i < 4; i++ ){
            floatType xi_i  = localPoints[ i ][ 0 ];
            floatType eta_i = localPoints[ i ][ 1 ];
            shapeFunctions[ i ] = 0.25 * ( 1 + xi_i * xi ) * xi_i * xi * ( 1 + eta_i * eta ) * eta_i * eta;
        }

        // Midpoint nodes
        shapeFunctions[ 4 ] = 0.5 * ( 1 -   xi * xi ) * ( 1 + localPoints[ 4 ][ 1 ] * eta ) * localPoints[ 4 ][ 1 ] * eta;
        shapeFunctions[ 5 ] = 0.5 * ( 1 - eta * eta ) * ( 1 + localPoints[ 5 ][ 0 ] *  xi ) * localPoints[ 5 ][ 0 ] *  xi;
        shapeFunctions[ 6 ] = 0.5 * ( 1 -   xi * xi ) * ( 1 + localPoints[ 6 ][ 1 ] * eta ) * localPoints[ 6 ][ 1 ] * eta;
        shapeFunctions[ 7 ] = 0.5 * ( 1 - eta * eta ) * ( 1 + localPoints[ 7 ][ 0 ] *  xi ) * localPoints[ 7 ][ 0 ] *  xi;

        // Center node
        shapeFunctions[ 8 ] = ( 1 - xi * xi) * ( 1 - eta * eta );

        return NULL;

    }

    errorOut evaluateGradQuadraticShapeFunctions( const floatType &xi, const floatType &eta, floatMatrix &gradShapeFunctions ){
        /*!
         * Evaluate the gradient of the shape functions w.r.t. the local coordinates
         * 
         * \param &xi: Local xi coordinate
         * \param &eta: Local eta coordinate
         * \param &gradShapeFunctions: The gradient of the shape functions w.r.t. the local coordinates
         */

        floatMatrix localPoints = {
                                      { -1, -1 },
                                      {  1, -1 },
                                      {  1,  1 },
                                      { -1,  1 },
                                      {  0, -1 },
                                      {  1,  0 },
                                      {  0,  1 },
                                      { -1,  0 },
                                      {  0,  0 }
                                  };

        gradShapeFunctions = floatMatrix( 9, floatVector( 2, 0 ) );

        // Corner nodes
        for ( unsigned int i = 0; i < 4; i++ ){
            floatType xi_i  = localPoints[ i ][ 0 ];
            floatType eta_i = localPoints[ i ][ 1 ];
            gradShapeFunctions[ i ][ 0 ] = 0.25 * xi_i * xi_i * xi * (1 + eta_i * eta) * eta_i * eta
                                         + 0.25 * (1 + xi_i * xi) * xi_i * (1 + eta_i * eta) * eta_i * eta;
            gradShapeFunctions[ i ][ 1 ] = 0.25 * (1 + xi_i * xi) * xi_i * xi * eta_i * eta_i * eta
                                         + 0.25 * (1 + xi_i * xi) * xi_i * xi * (1 + eta_i * eta) * eta_i;
        }

        // Midpoint nodes
        gradShapeFunctions[ 4 ][ 0 ] = -xi * (1 + localPoints[4][1] * eta) * localPoints[4][1] * eta;
        gradShapeFunctions[ 5 ][ 0 ] = 0.5 * (1 - eta * eta) * localPoints[5][0] * localPoints[5][0] *  xi
                                     + 0.5 * (1 - eta * eta) * (1 + localPoints[5][0] *  xi) * localPoints[5][0];
        gradShapeFunctions[ 6 ][ 0 ] = -xi * (1 + localPoints[6][1] * eta) * localPoints[6][1] * eta;
        gradShapeFunctions[ 7 ][ 0 ] = 0.5 * (1 - eta * eta) * localPoints[7][0] * localPoints[7][0] *  xi
                                     + 0.5 * (1 - eta * eta) * (1 + localPoints[7][0] *  xi) * localPoints[7][0];
        gradShapeFunctions[ 4 ][ 1 ] = 0.5 * (1 -   xi * xi) * localPoints[4][1] * localPoints[4][1] * eta
                                     + 0.5 * (1 -   xi * xi) * (1 + localPoints[4][1] * eta) * localPoints[4][1];
        gradShapeFunctions[ 5 ][ 1 ] = -eta * (1 + localPoints[5][0] *  xi) * localPoints[5][0] *  xi;
        gradShapeFunctions[ 6 ][ 1 ] = 0.5 * (1 -   xi * xi) * localPoints[6][1] * localPoints[6][1] * eta
                                     + 0.5 * (1 -   xi * xi) * (1 + localPoints[6][1] * eta) * localPoints[6][1];
        gradShapeFunctions[ 7 ][ 1 ] = -eta * (1 + localPoints[7][0] *  xi) * localPoints[7][0] *  xi;

        // Center node
        gradShapeFunctions[ 8 ][ 0 ] = -2 * xi * (1 - eta * eta);
        gradShapeFunctions[ 8 ][ 1 ] = -2 * (1 - xi * xi) * eta;

        return NULL;

    }

    errorOut interpolateFunction( const floatType &xi, const floatType &eta, floatMatrix &nodalValues, floatVector &answer ){
        /*!
         * Interpolate a function using the quadratic shape functions
         * 
         * \param &xi: Local xi coordinate
         * \param &eta: Local eta coordinate
         * \param &nodalValues: The values of the function at the nodes
         * \param &answer: The interpolated function
         */

        floatVector Ns;

        errorOut error;

        error = evaluateQuadraticShapeFunctions( xi, eta, Ns );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the shape functions" );
            result->addNext( error );
            return result;

        }

        if ( nodalValues.size( ) < 1 ){

            return new errorNode( __func__, "The nodal values have no entries." );

        }

        if ( nodalValues[ 0 ].size( ) != Ns.size( ) ){

            std::string message = "The nodal values must be of shape N outputs x N shape functions. The shape is:\n  nodalValues: " + std::to_string( nodalValues.size( ) ) + " x " + std::to_string( nodalValues[ 0 ].size( ) ) + "\n  shapeFunctions.size( ): " + std::to_string( Ns.size( ) );
            return new errorNode( __func__, message );

        }

        answer = vectorTools::dot( nodalValues, Ns );

        return NULL;

    }

    errorOut localGradientFunction( const floatType &xi, const floatType &eta, floatMatrix &nodalValues, floatMatrix &answer ){
        /*!
         * Compute the gradient of the interpolated function w.r.t. the local coordinates
         * 
         * \param &xi: Local xi coordinate
         * \param &eta: Local eta coordinate
         * \param &nodalValues: The values of the function at the nodes
         * \param &answer: The resulting gradient of the function w.r.t. the local coordinates
         */

        floatMatrix dNdxi;

        errorOut error;

        error = evaluateGradQuadraticShapeFunctions( xi, eta, dNdxi );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the shape functions" );
            result->addNext( error );
            return result;

        }

        if ( nodalValues.size( ) < 1 ){

            return new errorNode( __func__, "The nodal values have no entries." );

        }

        if ( nodalValues[ 0 ].size( ) != dNdxi.size( ) ){

            std::string message = "The nodal values must be of shape N outputs x N shape functions. The shape is:\n  nodalValues: " + std::to_string( nodalValues.size( ) ) + " x " + std::to_string( nodalValues[ 0 ].size( ) ) + "\n  shapeFunctions.size( ): " + std::to_string( dNdxi.size( ) );
            return new errorNode( __func__, message );

        }

        answer = vectorTools::dot( nodalValues, dNdxi );

        return NULL;

    }

    errorOut localJacobian( const floatType &xi, const floatType &eta, floatMatrix &nodalPositions, floatType &jacobian ){
        /*!
         * Compute the local jacobian of the quadratic element
         * 
         * \param &xi: Local xi coordinate
         * \param &eta: Local eta coordinate
         * \param &nodalPositions: The nodal positions
         * \param &jacobian: The surface jacobian
         */

        floatMatrix dxdxi;

        errorOut error;

        error = localGradientFunction( xi, eta, nodalPositions, dxdxi );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the local gradient" );
            result->addNext( error );
            return result;

        }

        if ( dxdxi.size( ) != 3 ){

            return new errorNode( __func__, "The nodal positions should be in 3d but are of dimension " + std::to_string( dxdxi.size( ) ) );

        }

        floatVector cross_product = { dxdxi[ 1 ][ 0 ] * dxdxi[ 2 ][ 1 ] - dxdxi[ 2 ][ 0 ] * dxdxi[ 1 ][ 1 ],
                                      dxdxi[ 2 ][ 0 ] * dxdxi[ 0 ][ 1 ] - dxdxi[ 0 ][ 0 ] * dxdxi[ 2 ][ 1 ],
                                      dxdxi[ 0 ][ 0 ] * dxdxi[ 1 ][ 1 ] - dxdxi[ 1 ][ 0 ] * dxdxi[ 0 ][ 1 ] };

        jacobian = vectorTools::l2norm( cross_product );

        return NULL;

    }

}
