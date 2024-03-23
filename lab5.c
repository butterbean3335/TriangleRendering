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
    char *output_filename = malloc(strlen(input_filename) + 5); // 4 for ".ppm", 1 for null terminator
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
    FILE *input_file = fopen(argv[1], "r");
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

    // Step 3: Allocate space for vertices and faces
    // (To be implemented)

    // Step 4: Read vertices and faces from input file
    // (To be implemented)

    // Step 5: Calculate bounding box
    // (To be implemented)

    // Step 6: Calculate E
    // (To be implemented)

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

    return 0;
}
