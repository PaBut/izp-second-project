/**Pavlo Butenko
 * xbuten00
 * Simple cluster analysis: 2D nearest neighbour.
 * Single linkage
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h> // sqrtf
#include <limits.h> // INT_MAX
#include <string.h>
#include <stdbool.h>

/*****************************************************************
 * Debugging macros. You can turn off their effect by defining the macro
 * NDEBUG, e.g.:
 * a) while compiling with the argument -DNDEBUG
 * b) in the file (on the line before #include <assert.h>
 *      #define NDEBUG
 */
#ifdef NDEBUG
#define debug(s)
#define dfmt(s, ...)
#define dint(i)
#define dfloat(f)
#else

// prints debugging string
#define debug(s) printf("- %s\n", s)

// prints formatted debugging output - used like printf
#define dfmt(s, ...) printf(" - "__FILE__":%u: "s"\n",__LINE__,__VA_ARGS__)

// prints variable debugging information - usage dint(variable)
#define dint(i) printf(" - " __FILE__ ":%u: " #i " = %d\n", __LINE__, i)

// prints debugging information about the variable with the type of float
// dfloat(variable)
#define dfloat(f) printf(" - " __FILE__ ":%u: " #f " = %g\n", __LINE__, f)

#endif

#define INT_MAX_DIGITS_COUNT 10

/*****************************************************************
 * Declaration of necessary data types:
 *
 * struct obj_t - object structure: identifier and coordinates
 * struct cluster_t - object cluster:
 * number of objects in the cluster,
 * cluster capacity (number of objects reserved for the cluster
 * in the array),
 * pointer to the cluster array.
 */

struct obj_t {
    int id;
    float x;
    float y;
};

struct cluster_t {
    int size;
    int capacity;
    struct obj_t *obj;
};

/*****************************************************************
 * Declaration of necessary functions.
 */

int return_error(char* msg, int return_value);

/*
 Initialization of cluster 'c'. Allocates memory for the object cap (capacity).
 The NULL pointer to the object array indicates a capacity of 0.
*/
int init_cluster(struct cluster_t *c, int cap)
{
    assert(c != NULL);
    assert(cap >= 0);

    if((c->obj = malloc(cap * sizeof(struct obj_t))) == NULL){
        return return_error("Allocation error", -1);
    }
    c->capacity = cap;
    c->size = 0;
    return 1;
}

/*
 Removes all cluster objects and initializes to an empty cluster.
 */
void clear_cluster(struct cluster_t *c)
{
    free(c->obj);
    c->size = 0;
    c->capacity = 0;
}


/*
 Changes the capacity of cluster 'c' to the 'new_cap'.
 */
struct cluster_t *resize_cluster(struct cluster_t *c, int new_cap)
{
    assert(c);
    assert(c->capacity >= 0);
    assert(new_cap >= 0);

    if (c->capacity >= new_cap)
        return c;

    size_t size = sizeof(struct obj_t) * new_cap;

    void *arr = realloc(c->obj, size);
    if (arr == NULL)
        return NULL;

    c->obj = (struct obj_t*)arr;
    c->capacity = new_cap;
    return c;
}

/*
 Adds object 'obj' to the end of the 'c' cluster.
 Expands the cluster if the object won't fit in it.
 */
int append_cluster(struct cluster_t *c, struct obj_t obj)
{
    if(c->size + 1 > c->capacity){
        if(resize_cluster(c, c->size+1) == NULL){
            return return_error("Allocation error", -1);
        }
    }
    c->obj[c->size++] = obj;
    return 1;
}

/*
 Sorts the objects in cluster 'c' in ascending order by their id.
 */
void sort_cluster(struct cluster_t *c);

/*
 Adds the 'c2' objects to the 'c1' cluster. The 'c1' cluster will be expanded if necessary.
 The objects in the 'c1' cluster will be sorted in ascending order by the id.
 The cluster 'c2' won't be changed.
 */
int merge_clusters(struct cluster_t *c1, struct cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c2 != NULL);

    //Checking here if the merging can fit in allocated memory for objects in c1
    if(c1->size + c2->size > c1->capacity){
        if(resize_cluster(c1, c1->size + c2->size) == NULL){
            return return_error("Allocation error", -1);
        }
    }
    for(int i = 0; i < c2->size; i++){
        c1->obj[c1->size + i] = c2->obj[i];
    }
    c1->size+=c2->size;
    sort_cluster(c1);
    return 1;
}

/**********************************************************************/
/* Working with a cluster array */

/*
 Removes a cluster from the 'carr' cluster array. The cluster array contains 'narr' elements
 (cluster). The cluster to remove is located at index 'idx'. The function returns a new
 number of clusters in the array.
*/
int remove_cluster(struct cluster_t *carr, int narr, int idx)
{
    assert(idx < narr);
    assert(narr > 0);

    clear_cluster(&carr[idx]);

    for(int i = idx; i < narr - 1; i++){
        carr[i] = carr[i+1];
    }
    return narr-1;
}

/*
Calculates the Euclidean distance between two objects.
 */
float obj_distance(struct obj_t *o1, struct obj_t *o2)
{
    assert(o1 != NULL);
    assert(o2 != NULL);

    return sqrt(pow(o1->x - o2->x, 2) + pow(o1->y - o2->y, 2));
}

/*
 It calculates the distance between two clusters.
*/
float cluster_distance(struct cluster_t *c1, struct cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c1->size > 0);
    assert(c2 != NULL);
    assert(c2->size > 0);

    float min_distance = obj_distance(&c1->obj[0], &c2->obj[0]);
    for(int i = 0; i < c1->size; i++){
        for(int j = 0; j < c2->size; j++){
            float tmp_distance = obj_distance(&c1->obj[i], &c2->obj[j]);
            if(tmp_distance < min_distance){
                min_distance = tmp_distance;
            }       
        }
    }
    return min_distance;
}

/*
 Finds the two closest clusters. In a cluster array 'carr' of size 'narr'
 searches for the two nearest clusters. It identifies the found clusters by their indixes in the array
 'carr'. The function stores the found clusters (indices in the 'carr' array) in memory in
 address 'c1' and 'c2' respectively.
*/
void find_neighbours(struct cluster_t *carr, int narr, int *c1, int *c2)
{
    assert(narr > 0);

    float min_distance = cluster_distance(&carr[0], &carr[1]);
    for(int i = 0; i < narr; i++){
        for(int j = i+1; j < narr; j++){
            int tmp_distance = cluster_distance(&carr[i], &carr[j]);
            if(tmp_distance <= min_distance){
                min_distance = tmp_distance;
                *c1 = i;
                *c2 = j;
            }
        }
    }
}

// helper function for cluster sorting
static int obj_sort_compar(const void *a, const void *b)
{
    const struct obj_t *o1 = (const struct obj_t *)a;
    const struct obj_t *o2 = (const struct obj_t *)b;
    if (o1->id < o2->id) return -1;
    if (o1->id > o2->id) return 1;
    return 0;
}

/*
 Sorts objects in the cluster in ascending order
  according to their identifiers.
*/
void sort_cluster(struct cluster_t *c)
{
    qsort(c->obj, c->size, sizeof(struct obj_t), &obj_sort_compar);
}

/*
 Prints cluster 'c' to stdout.
*/
void print_cluster(struct cluster_t *c)
{
    for (int i = 0; i < c->size; i++)
    {
        if (i) putchar(' ');
        printf("%d[%g,%g]", c->obj[i].id, c->obj[i].x, c->obj[i].y);
    }
    putchar('\n');
}

void free_clusters(struct cluster_t* arr, int len);
bool convert_int(char* str, int* result, bool must_be_higher_zero);
bool check_if_unique(struct cluster_t* arr, int len);
int return_error_with_file_closing(FILE* file, char* msg, int return_value);
int return_error_with_deallocation_file_closing(FILE* file, struct cluster_t* arr,
 int count, char* msg, int return_value);
int return_error_with_deallocation(struct cluster_t* arr,
 int count, char* msg, int return_value);
/*
 Retrieves objects from the file 'filename'. 
 For each object it creates a cluster and saves it into the cluster array.
 Allocates space for the array of all clusters and a pointer to the first
 position of the array (the pointer to the first cluster in the allocated array)
 is stored in memory, where the parameter 'arr' is referenced.
 The function returns the number of loaded objects (clusters).
 In case of an error, it stores a NULL value in the memory referenced by 'arr'.
*/
int load_clusters(char *filename, struct cluster_t **arr)
{
    //maximal count of digits in int plus one digit just to check if it's int plus '\0'
    const int STR_COUNT_ARG_LEN = 6 + 1 + INT_MAX_DIGITS_COUNT + 1 + 1; 

    assert(arr != NULL);

    FILE* file;
    file = fopen(filename, "r");
    if(file == NULL){
        return return_error("Argument error\n%s doesn\'t exist\n", 0);
    }
    char buffer[STR_COUNT_ARG_LEN];
    if(buffer == NULL){
        return return_error_with_file_closing(file, "Allocation error", 0);
    }

    fscanf(file, "%18s", buffer);

    //Checking here if a file in a proper format(must contain count= in the beginning)
    if(strcmp(strtok(buffer, "="), "count") != 0){
        return return_error_with_file_closing(file, "Argument error\nArgument file isn\'t in a proper format", 0);
    }
    strcpy(buffer, strtok(NULL, "="));

    int count;
    //Checking here if count argument in file is an integer higher than 0
    if(!convert_int(buffer, &count, true)){
        return return_error_with_file_closing(file, 
        "Argument file error\nThe value of [count] must be integer higher than 0", 0);
    }

    *arr = malloc(count * sizeof(struct cluster_t));
    if(arr == NULL){
        return return_error_with_file_closing(file, "Allocation error", 0);
    }

    int i = 0;
    //Setting file values to clusters
    for(; i < count; i++){
        int property_array[3];
        //Parsing values in the line and checking if they are int
        for(int j = 0; j < 3; j++){
            if(fscanf(file, "%18s", buffer) == EOF){
            return return_error_with_deallocation_file_closing(file, *arr, i, 
            "Argument file error\nCount of objects given \
in the beginning of the file is higher than initial count of objects in the file", 0);

            }

            if(!convert_int(buffer, &property_array[j], false)){
                return return_error_with_deallocation_file_closing(file, *arr, i, 
                "Argument file error\nProperies of objects isn\'t proper, they must be int type", 0);
            }
        }

        struct obj_t obj = {property_array[0], floor(property_array[1]), floor(property_array[2])};
        //Checking if there properties of multiple objects in a line 
        int last_symbol = getc(file);
        if(!(last_symbol == '\n' || last_symbol == '\r' || last_symbol == EOF)){
            return return_error_with_deallocation_file_closing(file, *arr, i, 
            "Argument file error\nProperies of more than one object can\'t be on a line", 0);
        }

        if(obj.x < 0 || obj.x >1000 || obj.y < 0 || obj.y > 1000){
            return return_error_with_deallocation_file_closing(file, *arr, i, 
            "Argument file error\nThe value of coordinates \
of some cluster must be in range between 0 and 1000 included", 0);

        }
        if(init_cluster(&(*arr)[i], 1) == -1){
            return return_error_with_deallocation_file_closing(file, *arr, i, "", 0);
        }
        if(append_cluster(&(*arr)[i], obj) == -1){
            return return_error_with_deallocation_file_closing(file, *arr, i, "", 0);
        }
        (*arr)[i].size = 1;
    }

    fclose(file);
    //Checking if all objects have unique id
    if(!check_if_unique(*arr, i)){
        return 0;
    }
    return i;
}

/*
 Tisk pole shluku. Parametr 'carr' je ukazatel na prvni polozku (shluk).
 Tiskne se prvnich 'narr' shluku.
*/
void print_clusters(struct cluster_t *carr, int narr)
{
    printf("Clusters:\n");
    for (int i = 0; i < narr; i++)
    {
        printf("cluster %d: ", i);
        print_cluster(&carr[i]);
    }
}


//Function to deallocate array of clusters
void free_clusters(struct cluster_t* arr, int len){
    for(int i = 0; i < len; i++){
        clear_cluster(&arr[i]);
    }
    free(arr);
}

bool check_if_unique(struct cluster_t* arr, int len){
    for(int j = 0; j < len; j++){
        for(int k = j + 1; k < len; k++){
            if(arr[j].obj[0].id == arr[k].obj[0].id){
                return return_error_with_deallocation(arr, len,
                 "Argument file error\nID in a cluster isn\'t unique", 0);
            }
        }
    }
    return true;
}

/*
Function to return error with message to stderr and 
deallocating the object with closing the file
*/
int return_error_with_deallocation_file_closing(FILE* file, struct cluster_t* arr,
 int count, char* msg, int return_value){
    free_clusters(arr, count);
    return return_error_with_file_closing(file, msg, return_value);
}

//Function to return error with message to stderr and deallocating the object
int return_error_with_deallocation(struct cluster_t* arr,
 int count, char* msg, int return_value){
    free_clusters(arr, count);
    return return_error(msg, return_value);
}

//Function to return error with message to stderr and closing the file
int return_error_with_file_closing(FILE* file, char* msg, int return_value){
    fclose(file);
    return return_error(msg, return_value);
}

//Function to return error with message to stderr
int return_error(char* msg, int return_value){
    fprintf(stderr, "%s\n", msg);
    return return_value;
}

//Function converts string to int and in the same time checks if it's an int
bool convert_int(char* str, int* result, bool must_be_higher_zero){
    int i = 0;
    unsigned max_digits_with_sign = INT_MAX_DIGITS_COUNT;
    //Number is negative, then maximal count of symbols in int number can be on one higher
    if(str[i] == '-'){
        i = 1;
        max_digits_with_sign++;
    }
    //Count of symbols in given string is higher than int number can store digits
    if(strlen(str) > max_digits_with_sign){
        return false;
    }
    //Checking if there're only digits in given string 
    for(; str[i] != '\0'; i++){
        if(str[i] < '0' || str[i] > '9'){
            return false;
        }
    }
    //Checking if number in int range
    long int temp = atol(str);
    if(temp > INT_MAX || temp < INT_MIN){
        return false;
    }
    if(must_be_higher_zero && temp <= 0){
        return false;
    }

    *result = (int)temp;
    return true;
}


int main(int argc, char *argv[])
{
    struct cluster_t *clusters;

    //Checking if arguments given to program are proper
    if(!(argc == 2 || argc == 3)){
        return return_error("Argument error\nWrong amount of arguments", EXIT_FAILURE);
    }

    int final_count = 1;

    if(argc == 3){
    	if(!convert_int(argv[2], &final_count, true)){
            return return_error("Argument error\nArgument [final_cluster_count] isn\'t in a proper format. \
It must be an integer higher than 0", EXIT_FAILURE);
        }
    }

    int len = load_clusters(argv[1], &clusters);

    //Checking if there was any error in load_clusters()
    if(len == 0){
        return EXIT_FAILURE;
    }
    else if(len < final_count){
        return return_error_with_deallocation(clusters, len,"Final count of the clusters \
can\'t be higher than count of clusters got from the input file", EXIT_FAILURE);
    }

    //nearest neighbors algorithm implemented in the cycle below
    while(len > final_count){
        int c1, c2;
        find_neighbours(clusters, len, &c1, &c2);
        if(merge_clusters(&clusters[c1], &clusters[c2]) == -1){
            return return_error_with_deallocation(clusters, len, "", EXIT_FAILURE);
        } 
        len = remove_cluster(clusters, len, c2);
    }
    print_clusters(clusters, len);

    free_clusters(clusters, len);

    return 0;
}
