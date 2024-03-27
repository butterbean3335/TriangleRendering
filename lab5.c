//triangle rendering lab

/*
 * 1) Setup
 *      -> check cmd line arguments, verify right number
 *      -> open inputfile "rb" for reading
 *      -> setup outputfilename using ppm extension
 *      -> open outputfile "wb" for writing
 * 2) Parse input file header
 *      -> number of verticies and number of faces present in input file, please
 *          print to ensure right numbers recieved
 * 3) Allocate space for vertices and faces
 * 4) You have space for verticies and faces, and know how many there are of each, read
 *      further in input file up to number of verticies or faces and store in allocated space
 * 5) calculate bounding box, minimum XYZ vector, max XYZ vector, center of these two vectors
 * 6) Calculate E, which is the largest difference between min and max in X, Y, or Z direction
 * 7) calculate camera position and orientation, default position {1, 0, 0}, up is {0, 0, 1}
 * 8) move and scale the camera
 * 9) determine 3D coordintes bounding the image, being the left, top, right, bottom, and topleft
 * 10) for each pixel, (r,c)
 *          -> calculate vector coordinates
 *          -> find plane equation, ABC
 *          -> find n and d
 *          -> find intersection
 *          -> determine if we see the triangle
 * 11) write ppm image
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define ROWS 256
#define COLS 256

void rotate (float objVector[3], float rotX[3][3], float rotY[3][3], float rotZ[3][3], float result[3])
{
    int row,col;
    float tmp1[3] = {0,0,0};
    float tmp2[3] = {0,0,0};
    // ROTATE X
    for(col = 0; col < 3; col++) {
        for(row = 0; row < 3; row++) {
            tmp1[col] += objVector[row] * rotX[row][col];
        }    
    }
    // ROTATE Y
    for(col = 0; col < 3; col++) {
        for(row = 0; row < 3; row++) {
            tmp2[col] += tmp1[row] * rotY[row][col];
        }    
    }
    // ROTATE Z
    for(col = 0; col < 3; col++) {
        for(row = 0; row < 3; row++) {
            result[col] += tmp2[row] * rotZ[row][col];
        }    
    }

}
// Function to calculate the largest extent E
float calculateExtent(float *min, float *max) {
    float extent[3]; // Store the extents in X, Y, and Z directions
    float E; // Largest extent of the bounding box

    // Calculate the extents in X, Y, and Z directions
    for (int i = 0; i < 3; i++) {
        extent[i] = max[i] - min[i];
    }

    // Find the largest extent among X, Y, and Z
    E = extent[0];
    for (int i = 1; i < 3; i++) {
        if (extent[i] > E) {
            E = extent[i];
        }
    }

    return E;
}

// Function to parse input file header
void parseHeader(FILE *input_file, int *num_vertices, int *num_faces) {
    char line[256];
    while (fgets(line, sizeof(line), input_file) != NULL) {
        if (strstr(line, "element vertex") != NULL) {
            sscanf(line, "element vertex %d", num_vertices);
        } else if (strstr(line, "element face") != NULL) {
            sscanf(line, "element face %d", num_faces);
        } else if (strcmp(line, "end_header\n") == 0) {
            break;  // End of header
        }
    }
}

// Function to generate output filename from input filename
char *getOutputFilename(const char *input_filename) {
    char *output_filename = malloc(strlen(input_filename) + 5); //4 for ".ppm", 1 for null terminator
    if (output_filename == NULL) {
        printf("Memory allocation error");
        exit(1);
    }
    strcpy(output_filename, input_filename);
    char *extension = strrchr(output_filename, '.');
    if (extension == NULL) {
        printf("Invalid input file name\n");
        free(output_filename);
        exit(1);
    }
    strcpy(extension, ".ppm");
    return output_filename;
}

// Function to calculate bounding box
void calculateBoundingBox(double *vertices, int num_vertices, float *min, float *max, float *center) {
    // Initialize min and max to the first vertex
    for (int i = 0; i < 3; i++) {
        min[i] = vertices[i * 3]; // Initialize with the value of the first vertex
        max[i] = vertices[i * 3]; // Initialize with the value of the first vertex
    }

    // Iterate through all vertices to find min and max along each axis
    for (int i = 0; i < num_vertices; i++) {
        // X axis
        if (vertices[i * 3] < min[0]) {
            min[0] = vertices[i * 3];
        }
        if (vertices[i * 3] > max[0]) {
            max[0] = vertices[i * 3];
        }
        // Y axis
        if (vertices[i * 3 + 1] < min[1]) {
            min[1] = vertices[i * 3 + 1];
        }
        if (vertices[i * 3 + 1] > max[1]) {
            max[1] = vertices[i * 3 + 1];
        }
        // Z axis
        if (vertices[i * 3 + 2] < min[2]) {
            min[2] = vertices[i * 3 + 2];
        }
        if (vertices[i * 3 + 2] > max[2]) {
            max[2] = vertices[i * 3 + 2];
        }
    }

    // Calculate center
    for (int i = 0; i < 3; i++) {
        center[i] = (min[i] + max[i]) / 2.0;
    }
}

void calculateCamera(float xAngle, float yAngle, float zAngle, float *camera, float *up) {
    float xRad = xAngle * (float)(M_PI/180);
    float yRad = yAngle * (float)(M_PI/180);
    float zRad = zAngle * (float)(M_PI/180);
    
    // Setup rotation matrices
    float rotateXMatrix[3][3] = {{1, 0, 0}, {0, cos(xRad), -sin(xRad)}, {0, sin(xRad), cos(xRad)}};
    float rotateYMatrix[3][3] = {{cos(yRad), 0, sin(yRad)}, {0, 1, 0}, {-sin(yRad), 0, cos(yRad)}};
    float rotateZMatrix[3][3] = {{cos(zRad), -sin(zRad), 0}, {sin(zRad), cos(zRad), 0}, {0, 0, 1}};


    float cameraRotated[3] = {0, 0, 0};
    float upRotated[3] = {0, 0, 0};
    int row, col;
    float tmp1[3] = {0,0,0};
    float tmp2[3] = {0,0,0};
    
    // Rotate camera vector
    for (col = 0; col < 3; col++) {
        for (row = 0; row < 3; row++) {
            tmp1[col] += camera[row] * rotateXMatrix[row][col];
        }    
    }
    for (col = 0; col < 3; col++) {
        for (row = 0; row < 3; row++) {
            tmp2[col] += tmp1[row] * rotateYMatrix[row][col];
        }    
    }
    for (col = 0; col < 3; col++) {
        for (row = 0; row < 3; row++) {
            cameraRotated[col] += tmp2[row] * rotateZMatrix[row][col];
        }    
    }
    
    rotate(up,rotateXMatrix,rotateYMatrix,rotateZMatrix,upRotated);
    up[0] = upRotated[0];
    up[1] = upRotated[1];
    up[2] = upRotated[2];
    
    // Update camera and up vectors
    for (int i = 0; i < 3; i++) {
        camera[i] = cameraRotated[i];
        up[i] = upRotated[i];
    }
}

float *sumVectors(float vector1[3], float vector2[3], float result[3]){
    for(int i=0; i<3; i++){
        result[i] = vector1[i]+vector2[i];
    }
}


void subVectors(float *vector1, float *vector2, float *result) {
    for (int i = 0; i < 3; i++) {
        result[i] = vector1[i] - vector2[i];
    }
}


float *multVector(float f, float vector[3]){
    static float result[3] = {0};
    for(int i=0; i<3; i++){
        result[i] = f * vector[i];
    }
    return result;
}

float *crossProduct(float vector1[3], float vector2[3], float result[3]){
    result[0] = (vector1[1] * vector2[2]) - (vector1[2]*vector2[1]);
    result[1] = (vector1[2] * vector2[0]) - (vector1[0]*vector2[2]);
    result[2] = (vector1[0] * vector2[1]) - (vector1[1]*vector2[0]);
}

float dotProduct(float vector1[3], float vector2[3]){
    float result;
    result = 0;
    result = ((vector1[0] * vector2[0]) + (vector1[1]*vector2[1]) + (vector2[1] * vector2[2]));
    return result;
}

float magnitude(float vec[3]) {
    return sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
}

// Function to determine 3D coordinates bounding the image
void determineImageBounds(float *center, float *camera, float *up, float E, float aspect_ratio, float *left, float *right, float *top, float *bottom, float *topleft) {
    // Calculate left and right vectors
    float temp[3];
    subVectors(center, camera, temp);
    crossProduct(up, temp, left);
    float a = magnitude(left);
    for (int i = 0; i < 3; i++) {
        left[i] = (E / (2 * a)) * left[i] + center[i];
    }
    crossProduct(temp, up, right);
    for (int i = 0; i < 3; i++) {
        right[i] = (E / (2 * a)) * right[i] + center[i];
    }

    // Calculate top and bottom vectors
    for (int i = 0; i < 3; i++) {
        top[i] = (E / 2) * up[i] + center[i];
        bottom[i] = (-1 * E / 2) * up[i] + center[i];
    }

    // Calculate top left corner
    for (int i = 0; i < 3; i++) {
        topleft[i] = (E / 2) * up[i] + left[i];
    }
}

void posZeros(float matrix[3][3]) 
{
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            if(matrix[i][j] > 0 && matrix[i][j] < 0.0001)
                matrix[i][j] = 0;
            if(matrix[i][j] < 0 && (0 - matrix[i][j]) < 0.0001 )
                matrix[i][j] *= -1;
            if(matrix[i][j] == 0)
                matrix[i][j] = 0;
        }
    }
}


int main(int argc, char *argv[]) {
    // Step 1: Setup
    // Check command-line arguments
    if (argc != 5) {
        printf("Usage: %s <input_file> X Y Z\n", argv[0]);
        return 1;
    }

    // Extract filename prefix from input file
    char *output_filename = getOutputFilename(argv[1]);

    // Open input file for reading
    FILE *input_file = fopen(argv[1], "rb");
    if (input_file == NULL) {
        printf("Error opening input file");
        free(output_filename);
        return 1;
    }

    // Open output file for writing
    FILE *output_file = fopen(output_filename, "wb");
    if (output_file == NULL) {
        printf("Error opening output file");
        fclose(input_file);
        free(output_filename);
        return 1;
    }

    // Step 2: Parse input file header
    int num_vertices = 0, num_faces = 0;
    parseHeader(input_file, &num_vertices, &num_faces);
    //printf("Number of vertices: %d\n", num_vertices);
    //printf("Number of faces: %d\n", num_faces);

    // Consume the newline character after parsing the header
    // Move the file pointer back to the current position
    fseek(input_file, 0, SEEK_SET);
    // Move the file pointer to the 10th line, standard for ply files
    int line_number = 0;
    char line[256];

    while (fgets(line, sizeof(line), input_file) != NULL) {
        line_number++;
        if (line_number == 9) {
            break;
        }
    }

    // Read and print characters from the current position until a newline or EOF is encountered
    /*char c;
    while ((c = fgetc(input_file)) != '\n' && c != EOF) {
        putchar(c);
    }
    putchar('\n'); 
    */

    // Step 3: Allocate space for vertices and faces
    double *vertices = (double *)calloc(num_vertices * 3, sizeof(double)); //stores __ __ __
    int *faces = (int *)calloc(num_faces * 4, sizeof(int));               //stores 3 __ __ __

    // Step 4: Read vertices and faces from input file
    // Read vertices
    int i;
    for (i = 0; i < num_vertices; i++) {
        fscanf(input_file, "%lf %lf %lf", &vertices[i * 3], &vertices[i * 3 + 1], &vertices[i * 3 + 2]);
        //printf("%lf %lf %lf\n", vertices[i * 3], vertices[i * 3 + 1], vertices[i * 3 + 2]);
    }

    // Read faces
    for (i = 0; i < num_faces; i++) {
        fscanf(input_file, "%d %d %d %d", &faces[i * 4], &faces[i * 4 + 1], &faces[i * 4 + 2], &faces[i * 4 + 3]);
        //printf("%d %d %d %d\n", faces[i * 4], faces[i * 4 + 1], faces[i * 4 + 2], faces[i * 4 + 3]);
    }



    // Step 5: Calculate bounding box
    // will have max, min, and center of verticies read into verticies array
    // need to store X, Y, and Z, and so min and max will be an array with size 3
    float min[3], max[3], center[3];
    calculateBoundingBox(vertices, num_vertices, min, max, center);

    // Print or use min, max, and center as needed
    printf("Bounding Box:\n");
    printf("Minimum X: %f\n", min[0]);
    printf("Minimum Y: %f\n", min[1]);
    printf("Minimum Z: %f\n", min[2]);
    printf("Maximum X: %f\n", max[0]);
    printf("Maximum Y: %f\n", max[1]);
    printf("Maximum Z: %f\n", max[2]);
    printf("Center of bounding box: (%f, %f, %f)\n", center[0], center[1], center[2]);


    // Step 6: Calculate E
    float E = calculateExtent(min, max);
    printf("Largest extent E: %f\n", E);


    // Step 7: Calculate camera position and orientation
    float camera[3] = {1, 0, 0};
    float up[3] = {0, 0, 1};
    calculateCamera(atof(argv[2]), atof(argv[3]), atof(argv[4]), camera, up);
    printf("Camera-> X: %f\tY: %f\tZ: %f\n",camera[0],camera[1],camera[2]);

    // Step 8: Move and scale the camera
    for(i = 0; i < 3; i++) {
        camera[i] = 1.5*E*camera[i] + center[i];
    }
    printf("Camera(move and scale)-> X: %f\tY: %f\tZ: %f\n",camera[0],camera[1],camera[2]);
    printf("Up(move and scale)-> X: %f\tY: %f\tZ: %f\n",up[0],up[1],up[2]);

    // Step 9: Determine 3D coordinates bounding the image
    float left[3], right[3], top[3], bottom[3], topleft[3];

    determineImageBounds(center, camera, up, E, ROWS/COLS, left, right, top, bottom, topleft);
    printf("camera left/right: \t\t(%0.1f %0.1f %0.1f) \t(%0.1f %0.1f %0.1f)\n", left[0], left[1], left[2], right[0], right[1], right[2]);
    printf("camera top/bottom/topleft: \t(%0.1f %0.1f %0.1f) \t\t(%0.1f %0.1f %0.1f) \t(%0.1f %0.1f %0.1f)\n", top[0], top[1], top[2], bottom[0], bottom[1], bottom[2], topleft[0], topleft[1], topleft[2]);
  
  
    // Step 10: Render pixels
    float top2[2], bottom2[2];
    top2[0] = top[0];
    top2[1] = top[1];
    top2[2] = top[2];
    bottom2[0] = bottom[0];
    bottom2[1] = bottom[1];
    bottom2[2] = bottom[2];
    int r,c;
    unsigned char *pixels;
    float A, B, C, D, n, d, dot1, dot2, dot3, cdiv, rdiv;
    float zbuffer, image[3], v0[3], v1[3], v2[3], intersect[3], cross1[3], cross2[3];
    pixels = (unsigned char *)calloc(ROWS*COLS,1);
                
    float diff1[3], diff2[3];
    int color, maxColor;
    printf("Rendering...\n");
    for(int p = 0; p < ROWS*COLS; p++) {
            r = p/ROWS;
            c = p%COLS;
            zbuffer = 9999999;
            if(c == 0) printf("%d ", r);
            for(i = 0; i < 3; i++) {
                image[i] = topleft[i] + ((float)c/(COLS-1))*(right[i]-left[i]) + ((float)r/(ROWS-1))*(bottom[i]-top[i]);
            }
            for(i = 0; i < num_faces; i++){
                //get triangle points
                v0[0] = vertices[faces[i*4+1]*3];
                v0[1] = vertices[faces[i*4+1]*3+1];
                v0[2] = vertices[faces[i*4+1]*3+2];
                v1[0] = vertices[faces[i*4+2]*3];
                v1[1] = vertices[faces[i*4+2]*3+1];
                v1[2] = vertices[faces[i*4+2]*3+2];
                v2[0] = vertices[faces[i*4+3]*3];
                v2[1] = vertices[faces[i*4+3]*3+1];
                v2[2] = vertices[faces[i*4+3]*3+2];
                //printf("v0: <%f, %f, %f>\n",v0[0],v0[1],v0[2]);
                //printf("v1: <%f, %f, %f>\n",v1[0],v1[1],v1[2]);
                //printf("v2: <%f, %f, %f> %d\n",v2[0],v2[1],v2[2],i);

                //find the plane equation
                A = (v1[1] - v0[1])*(v2[2] - v0[2]) - (v1[2] - v0[2])*(v2[1] - v0[1]);
                B = (v1[2] - v0[2])*(v2[0] - v0[0]) - (v1[0] - v0[0])*(v2[2] - v0[2]);
                C = (v1[0] - v0[0])*(v2[1] - v0[1]) - (v1[1] - v0[1])*(v2[0] - v0[0]);
                D = -1*A*v0[0] + -1*B*v0[1] + -1*C*v0[2];

                n = -1*A*camera[0] + -1*B*camera[1] + -1*C*camera[2] - D;
                d = A*(image[0] - camera[0]) + B*(image[1] - camera[1]) + C*(image[2] - camera[2]);
                if(fabs(d) < 0.01) {
                    continue; 
                }

                intersect[0] = camera[0] + (n/d)*(image[0] - camera[0]);
                intersect[1] = camera[1] + (n/d)*(image[1] - camera[1]);
                intersect[2] = camera[2] + (n/d)*(image[2] - camera[2]);

                subVectors(v2, v0, diff1);
                subVectors(v1, v0, diff2);
                crossProduct(diff1, diff2, cross1);
                subVectors(intersect, v0, diff1);
                subVectors(v1, v0, diff2);
                crossProduct(diff1, diff2, cross2);
                dot1 = cross1[0]*cross2[0] + cross1[1]*cross2[1] + cross1[2]*cross2[2];
                //if(r == 0 && c == 0) printf("dot1: %f\n", dot1);

                subVectors(v0, v1, diff1);
                subVectors(v2, v1, diff2);
                crossProduct(diff1, diff2, cross1);
                subVectors(intersect, v1, diff1);
                crossProduct(diff1, diff2, cross2);
                dot2 = cross1[0]*cross2[0] + cross1[1]*cross2[1] + cross1[2]*cross2[2];
                //if(r == 0 && c == 0) printf("dot2: %f\n", dot2);

                subVectors(v1, v2, diff1);
                subVectors(v0, v2, diff2);
                crossProduct(diff1, diff2, cross1);
                subVectors(intersect, v2, diff1);
                crossProduct(diff1, diff2, cross2);
                dot3 = cross1[0]*cross2[0] + cross1[1]*cross2[1] + cross1[2]*cross2[2];
                //if(r == 0 && c == 0) printf("dot3: %f\n", dot3);
                
                //if(r == 0 && c == 0) printf("DOT 1: %f, DOT 2: %f, DOT 3: %f\n", dot1, dot2, dot3);

                if(dot1 < 0 || dot2 < 0 || dot3 < 0) {
                    continue; // skip triangle
                } 
                
                if( zbuffer < (n/d))
                    zbuffer = (n/d);
                pixels[r*COLS+c] = 155 + (i%100);
            }
        }



    // Step 11: Write ppm image
    fprintf(output_file,"P5 %d %d 255\n",COLS,ROWS);
    fwrite(pixels,COLS*ROWS,sizeof(unsigned char),output_file);


    // Close files
    fclose(input_file);
    free(output_filename);
    free(vertices);
    free(faces);
    fclose(output_file);
    return 0;
}
