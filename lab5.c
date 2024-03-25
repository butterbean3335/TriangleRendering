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

int main(int argc, char *argv[]) {
    // Step 1: Setup
    // Check command-line arguments
    if (argc != 2) {
        printf("Usage: %s <input_file>\n", argv[0]);
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
    printf("Number of vertices: %d\n", num_vertices);
    printf("Number of faces: %d\n", num_faces);

    // Consume the newline character after parsing the header
    // Move the file pointer back to the current position
    fseek(input_file, 0, SEEK_SET);
    // Move the file pointer to the 10th line
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
    putchar('\n'); // Print a newline at the end for formatting
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

    // Consume the newline character after reading vertices

    // Read faces
    for (i = 0; i < num_faces; i++) {
        fscanf(input_file, "%d %d %d %d", &faces[i * 4], &faces[i * 4 + 1], &faces[i * 4 + 2], &faces[i * 4 + 3]);
        //printf("%d %d %d %d\n", faces[i * 4], faces[i * 4 + 1], faces[i * 4 + 2], faces[i * 4 + 3]);
    }



    // Step 5: Calculate bounding box
    // will have max, min, and center of verticies read into verticies array
    // need to store X, Y, and Z, and so min and max will be an array with size 3
    float min[3];
    float max[3];
    float center[3];
    //step through all vertices to find max in X, Y, Z independently
    //step through all vertices to find min in X,Y,Z independently
    //center splits the difference between min and max
    // Initialize min and max to the first vertex
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

    // Print or use min, max, and center as needed
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
    // (To be implemented)

    // Step 8: Move and scale the camera
    // (To be implemented)

    // Step 9: Determine 3D coordinates bounding the image
    // (To be implemented)

    // Step 10: Render pixels
    // (To be implemented)

    // Step 11: Write ppm image
    // (To be implemented)

    // Close files
    fclose(input_file);
    fclose(output_file);
    free(output_filename);
    free(vertices);
    free(faces);

    return 0;
}
