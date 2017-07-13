#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <png.h>

#define COLOR_BIT_DEPTH 7
#define CHANNEL_SIZE (1 << COLOR_BIT_DEPTH)
#define IMAGE_X 2048
#define IMAGE_Y 1024

#define N_TL 0
#define N_T 1
#define N_TR 2
#define N_L 3
#define N_R 4
#define N_BL 5
#define N_B 6
#define N_BR 7

#define KD_LEAF_COUNT 128
#define KD_RED 0
#define KD_GREEN 1
#define KD_BLUE 2

#define STARTING_SEEDS 1

struct pixel_t {
    unsigned char r;
    unsigned char g;
    unsigned char b;
    unsigned char a;
};
typedef struct pixel_t pixel;

struct pixel_idx_t {
    pixel p;
    int idx;
};
typedef struct pixel_idx_t pixel_idx;


struct kd_pixel_leaf_t {
    pixel_idx pixels[KD_LEAF_COUNT];
    unsigned char count;
};
typedef struct kd_pixel_leaf_t kd_pixel_leaf;

struct kd_pixel_node_t {
    struct kd_pixel_node_t* left;
    struct kd_pixel_node_t* right;
    kd_pixel_leaf* leaf;
    unsigned int split_axis;
    unsigned char pos[3];
};
typedef struct kd_pixel_node_t kd_pixel_node;

struct kd_pixel_tree_t {
    int pixel_count;
    int leaf_count;
    int node_count;
    kd_pixel_node* root;
};
typedef struct kd_pixel_tree_t kd_pixel_tree;

int color_diff_sq(pixel c1, pixel c2) {
    int r = c1.r - c2.r;
    int g = c1.g - c2.g;
    int b = c1.b - c2.b;
    return r * r + g * g + b * b;
}

kd_pixel_node* create_kd_node() {
    kd_pixel_node* node = malloc(sizeof(kd_pixel_node));
    memset(node, 0, sizeof(kd_pixel_node));
    return node;
}

void destroy_kd_node(kd_pixel_node* node) {
    free(node);
}

kd_pixel_leaf* create_kd_leaf() {
    kd_pixel_leaf* leaf = malloc(sizeof(kd_pixel_leaf));
    memset(leaf, 0, sizeof(kd_pixel_leaf));
    return leaf;
}

void destroy_kd_leaf(kd_pixel_leaf* leaf) {
    free(leaf);
}

kd_pixel_tree* create_kd_tree() {
    kd_pixel_tree* tree = malloc(sizeof(kd_pixel_tree));
    memset(tree, 0, sizeof(kd_pixel_tree));
    tree->root = create_kd_node();
    tree->root->leaf = create_kd_leaf();
    tree->leaf_count = 1;
    return tree;
}

void destroy_kd_pixel_tree_r(kd_pixel_node* root) {
    if (root->leaf) {
        destroy_kd_leaf(root->leaf);
    } else {
        destroy_kd_pixel_tree_r(root->left);
        destroy_kd_pixel_tree_r(root->right);
    }
    destroy_kd_node(root);
}

void destroy_kd_tree(kd_pixel_tree* tree) {
    destroy_kd_pixel_tree_r(tree->root);
    free(tree);
}

kd_pixel_node* split_kd_leaf(kd_pixel_node* leaf_node) {
    /* find splitting plane between pixels in leaf */
    kd_pixel_leaf* leaf = leaf_node->leaf;
    int r = 0;
    int g = 0;
    int b = 0;
    int r_diff = 0;
    int g_diff = 0;
    int b_diff = 0;
    int split = KD_RED;
    int split_val;
    int i;
    for (i = 0; i < leaf->count; i++) {
        r += leaf->pixels[i].p.r;
        g += leaf->pixels[i].p.g;
        b += leaf->pixels[i].p.b;
    }
    r = (r / leaf->count) + (r % leaf->count != 0);
    g = (g / leaf->count) + (g % leaf->count != 0);
    b = (b / leaf->count) + (b % leaf->count != 0);
    for (int i = 0; i < leaf->count; i++) {
        r_diff += r <= leaf->pixels[i].p.r ? 1 : -1;
        g_diff += g <= leaf->pixels[i].p.g ? 1 : -1;
        b_diff += b <= leaf->pixels[i].p.b ? 1 : -1;
    }
    if (r_diff < 0) r_diff *= -1;
    if (g_diff < 0) g_diff *= -1;
    if (b_diff < 0) b_diff *= -1;
    split = r_diff < g_diff ? KD_RED : KD_GREEN;
    int tmp = r_diff < g_diff ? r_diff : g_diff;
    split = tmp < b_diff ? split : KD_BLUE;

    int count = leaf->count;
    kd_pixel_leaf* leaf_l = leaf;
    kd_pixel_leaf* leaf_r = create_kd_leaf();
    leaf_l->count = 0;
    leaf_r->count = 0;
    for (i = 0; i < count; i++) {
        kd_pixel_leaf* l;
        if (split == KD_RED) l = leaf->pixels[i].p.r < r ? leaf_l : leaf_r;
        else if (split == KD_GREEN) l = leaf->pixels[i].p.g < g ? leaf_l : leaf_r;
        else if (split == KD_BLUE) l = leaf->pixels[i].p.b < b ? leaf_l : leaf_r;
        l->pixels[l->count++] = leaf->pixels[i];
    }
    if (split == KD_RED) split_val = r;
    else if (split == KD_GREEN) split_val = g;
    else if (split == KD_BLUE) split_val = b;
    kd_pixel_node* node_l = create_kd_node();
    kd_pixel_node* node_r = create_kd_node();
    node_l->leaf = leaf_l;
    node_r->leaf = leaf_r;
    kd_pixel_node* top_node = leaf_node;
    top_node->leaf = NULL;
    top_node->left = node_l;
    top_node->right = node_r;
    top_node->split_axis = split;
    top_node->pos[split] = split_val;
    return top_node;
}

kd_pixel_node* find_subspace(kd_pixel_node* root, const pixel p) {
    kd_pixel_node* node = root;
    while (node->leaf == NULL) {
        int query_val = p.r;
        int split_val = node->pos[node->split_axis];
        if (node->split_axis == KD_RED) query_val = p.r;
        else if (node->split_axis == KD_GREEN) query_val = p.g;
        else if (node->split_axis == KD_BLUE) query_val = p.b;
        node = query_val < split_val ? node->left : node->right;
    }
    return node;
}

void add_pixel(kd_pixel_tree* tree, const pixel_idx p) {
    kd_pixel_node* node = find_subspace(tree->root, p.p);
    kd_pixel_leaf* leaf = node->leaf;
    if (leaf->count >= KD_LEAF_COUNT) {
        split_kd_leaf(node);
        node = find_subspace(node, p.p);
        leaf = node->leaf;
        tree->leaf_count++;
        if (leaf->count >= KD_LEAF_COUNT) {
            printf("LEAF SPLIT DIDN'T SPLIT!\n");
            abort();
        }
    }
    leaf->pixels[leaf->count++] = p;
    tree->pixel_count++;
}

void remove_pixel(kd_pixel_tree* tree, const pixel_idx p) {
    kd_pixel_node* node = find_subspace(tree->root, p.p);
    kd_pixel_leaf* leaf = node->leaf;
    int i = 0;
    for (i = 0; i < leaf->count; i++) {
        pixel_idx q = leaf->pixels[i];
        if (p.idx == q.idx) {
            #ifdef DEBUG
            if (p.p.r != q.p.r || p.p.g != q.p.g || p.p.b != q.p.b) {
                static int t = 0;
                printf("WTF?\n");
                printf("Idx: %d\n", p.idx);
                printf("%d\n", ++t);
            }
            #endif
            memmove(leaf->pixels + i, leaf->pixels + i + 1, (leaf->count - i - 1) * sizeof(pixel_idx));
            leaf->count--;
            tree->pixel_count--;
            break;
        }
    }
}

int closest_in_kd_leaf(const kd_pixel_leaf* leaf, const pixel query_pixel, pixel_idx* out) {
    int i = 0;
    int d_min = INT_MAX;
    for (i = 0; i < leaf->count; i++) {
        pixel_idx p = leaf->pixels[i];
        int d = color_diff_sq(p.p, query_pixel);
        if (d < d_min) {
            d_min = d;
            *out = p;
        }
    }
    return d_min;
}

int search_kd_pixel_tree(const kd_pixel_node* node, const pixel p, pixel_idx* out) {
    if (node->leaf) {
        return closest_in_kd_leaf(node->leaf, p, out);
    } else {
        int d = p.r - node->pos[KD_RED];
        kd_pixel_node *first, *second;
        if (node->split_axis == KD_RED) d = p.r - node->pos[KD_RED];
        else if (node->split_axis == KD_GREEN) d = p.g - node->pos[KD_GREEN];
        else if (node->split_axis == KD_BLUE) d = p.b - node->pos[KD_BLUE];
        first = d < 0 ? node->left : node->right;
        second = d < 0 ? node->right : node->left;
        int dist = search_kd_pixel_tree(first, p, out);
        if (d < 0) d *= -1;
        if (dist <= d) {
            return dist;
        } else {
            pixel_idx alt_out;
            int dist2 = search_kd_pixel_tree(second, p, &alt_out);
            if (dist2 < dist) {
                *out = alt_out;
                return dist2;
            } else {
                return dist;
            }
        }
    }
}

int find_closest(const kd_pixel_tree* kdtree, const pixel p, pixel_idx* out) {
    return search_kd_pixel_tree(kdtree->root, p, out);
}

pixel int_to_pixel(int x) {
    pixel c = {(x >> 16) & 0xFF,
        (x >> 8) & 0xFF,
        x & 0xFF};
    return c;
}

int write_png(const pixel* data, const char* path,
        const int size_x, const int size_y) {
    int status = -1;
    FILE* f;
    png_structp png = NULL;
    png_infop info = NULL;
    png_byte** rows = NULL;
    int pixel_size = 3;
    int depth = 8;

    /* acquire resources */
    f = fopen(path, "wb");
    if (f == NULL) goto open_fail;
    png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png == NULL) goto png_struct_fail;
    info = png_create_info_struct(png);
    if (info == NULL) goto info_struct_fail;
    /* setjmp for error handling */
    if (setjmp(png_jmpbuf(png))) goto png_failure;

    png_set_IHDR(png, info, size_x, size_y, depth,
            PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    rows = png_malloc(png, size_y * sizeof(png_byte*));
    size_t y, x;
    const pixel* ptr_data = data;
    for (y = 0; y < size_y; y++) {
        png_byte *row = png_malloc(png, size_x * pixel_size);
        rows[y] = row;
        for (x = 0; x < size_x; x++) {
            *row++ = ptr_data->r;
            *row++ = ptr_data->g;
            *row++ = ptr_data->b;
            ptr_data++;
        }
    }
    png_init_io(png, f);
    png_set_rows(png, info, rows);
    png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);

    status = 0;
png_failure:
info_struct_fail:
    png_destroy_write_struct(&png, &info);
png_struct_fail:
    fclose(f);
open_fail:
    return status;
}

int rand_in_range(int min, int max) {
    return min + rand() / (RAND_MAX / (max - min + 1) + 1);
}

void shuffle(pixel* pixels, const int size) {
    int i, j;
    for (i = 0; i < size - 1; i++) {
        j = rand_in_range(i, size-1);
        pixel tmp = pixels[i];
        pixels[i] = pixels[j];
        pixels[j] = tmp;
    }
}

int hue_sort(const void* p1, const void* p2) {
    /* sort based on hue angle */
    pixel* x = (pixel*) p1;
    pixel* y = (pixel*) p2;
    int hxn = x->g - x->b;
    int hxd = 2.0 * x->r - x->g - x->b;
    int hyn = y->g - y->b;
    int hyd = 2.0 * y->r - y->g - y->b;
    if (hxn < 0 && hyn >= 0) return -1;
    if (hxn >= 0 && hyn < 0) return 1;
    if (hxn < 0 && hyn < 0) {
        if (hxd < 0 && hyd >= 0) return -1;
        if (hxd >=0 && hyd < 0) return 1;
    }
    if (hxn > 0 && hyn > 0) {
        if (hxd < 0 && hyd >= 0) return 1;
        if (hxd >=0 && hyd < 0) return -1;
    }

    /* sort out zeroes in numerator */
    if (hxn == 0 && hxd < 0) {
        return hyn == 0 ? 0 : 1;
    }
    if (hyn == 0 && hyd < 0) {
        return hxn == 0 ? 0 : -1;
    }
    if ((hxn == 0 && hxd >= 0) || (hyn == 0 && hyd >= 0)) {
        return (hxn < hyn) ? -1 : (hxn > hyn);
    }

    /*
     * atan(v) is strictly increasing for v in the same quadrant
     * so only need to compare hxn/hxd < hyn/hyd, after handling zero case
     * note n1/d1 < n2/d2 -> n1*d2 < n2*d2 so long as d1 and d2 have same sign
     */
    int lhs = hxn * hyd;
    int rhs = hyn * hxd;
    return (lhs < rhs) ? -1 : (lhs > rhs);

}

int hsp_sort(const void* p1, const void* p2) {
    /*
     * copied from alienryderflex.com/hsp.html
     * appears to be based on Photoshop RGB->Grayscale conversion
     */
    pixel* x = (pixel*) p1;
    pixel* y = (pixel*) p2;
    float bx = 0.299 * x->r * x->r + 0.587 * x->g * x->g + 0.144 * x->b * x->b;
    float by = 0.299 * y->r * y->r + 0.587 * y->g * y->g + 0.144 * y->b * y->b;
    return (bx < by) ? -1 : (bx > by);
}

void color_sort(pixel* pixels, const int size) {
    qsort(pixels, size, sizeof(*pixels), &hue_sort);
}

void get_neighbour_idxs(const int idx, const int size_x, const int size_y, int* others) {
    int l_edge = idx % size_x == 0;
    int r_edge = idx % size_x == size_x - 1;
    int t_edge = idx < size_x;
    int b_edge = idx >= (size_x * (size_y - 1));
    int invalid = -1;
    others[N_TL] = !l_edge && !t_edge ? idx - size_x - 1 : invalid;
    others[N_T] = !t_edge ? idx - size_x : invalid;
    others[N_TR] = !r_edge && !t_edge ? idx - size_x + 1 : invalid;
    others[N_L] = !l_edge ? idx - 1 : invalid;
    others[N_R] = !r_edge ? idx + 1 : invalid;
    others[N_BL] = !l_edge && !b_edge ? idx - 1 + size_x : invalid;
    others[N_B] = !b_edge ? idx + size_x : invalid;
    others[N_BR] = !r_edge && !b_edge ? idx + 1 + size_x : invalid;
}

int get_available_neighbour_idxs(const int idx, const pixel* data, const int size_x, const int size_y, int* others) {
    int l_edge = idx % size_x == 0;
    int r_edge = idx % size_x == size_x - 1;
    int t_edge = idx < size_x;
    int b_edge = idx >= (size_x * (size_y - 1));
    int invalid = -1;
    int* ptr = others;
    *ptr = !l_edge && !t_edge ? idx - size_x - 1 : invalid;
    ptr += (*ptr == invalid || data[*ptr].a == 0x1) ? 0 : 1;
    *ptr = !t_edge ? idx - size_x : invalid;
    ptr += (*ptr == invalid || data[*ptr].a == 0x1) ? 0 : 1;
    *ptr = !r_edge && !t_edge ? idx - size_x + 1 : invalid;
    ptr += (*ptr == invalid || data[*ptr].a == 0x1) ? 0 : 1;
    *ptr = !l_edge ? idx - 1 : invalid;
    ptr += (*ptr == invalid || data[*ptr].a == 0x1) ? 0 : 1;
    *ptr = !r_edge ? idx + 1 : invalid;
    ptr += (*ptr == invalid || data[*ptr].a == 0x1) ? 0 : 1;
    *ptr = !l_edge && !b_edge ? idx - 1 + size_x : invalid;
    ptr += (*ptr == invalid || data[*ptr].a == 0x1) ? 0 : 1;
    *ptr = !b_edge ? idx + size_x : invalid;
    ptr += (*ptr == invalid || data[*ptr].a == 0x1) ? 0 : 1;
    *ptr = !r_edge && !b_edge ? idx + 1 + size_x : invalid;
    ptr += (*ptr == invalid || data[*ptr].a == 0x1) ? 0 : 1;
    return ptr - others;
}

void get_neighbours(const int idx, const pixel* data,
        const int size_x, const int size_y,
        pixel* others) {

    int idxs[8];
    pixel invalid = {0, 0, 0, 0x7F};

    get_neighbour_idxs(idx, size_x, size_y, idxs);
    others[N_TL] = idxs[N_TL] > -1 ? data[idxs[N_TL]] : invalid;
    others[N_T] = idxs[N_T] > -1 ? data[idxs[N_T]] : invalid;
    others[N_TR] = idxs[N_TR] > -1 ? data[idxs[N_TR]] : invalid;
    others[N_L] = idxs[N_L] > -1 ? data[idxs[N_L]] : invalid;
    others[N_R] = idxs[N_R] > -1 ? data[idxs[N_R]] : invalid;
    others[N_BL] = idxs[N_BL] > -1 ? data[idxs[N_BL]] : invalid;
    others[N_B] = idxs[N_B] > -1 ? data[idxs[N_B]] : invalid;
    others[N_BR] = idxs[N_BR] > -1 ? data[idxs[N_BR]] : invalid;
}

kd_pixel_tree* build_new_tree(const pixel* buffer, const int size_x, const int size_y) {
    kd_pixel_tree* tree = create_kd_tree();
    int idx = -1;
    while (++idx < size_x * size_y) {
        pixel p = buffer[idx];
        int neighbours[8] = {0};
        if (p.a == 1 && get_available_neighbour_idxs(idx, buffer, size_x, size_y, neighbours)) {
            pixel_idx pdx = {p, idx};
            add_pixel(tree, pdx);
        }
    }
    return tree;
}

void fill(const pixel* pixels, pixel* buffer,
        const int size_x, const int size_y) {
    int i = 0;
    const int size = size_x * size_y;
    const pixel* ptr = pixels;
    kd_pixel_tree* kdtree = create_kd_tree();
    for (i = 0 ; i < STARTING_SEEDS; i++) {
        /* TODO: use different ranges for seeds */
        int idx = rand_in_range(0, size - 1);
        buffer[idx] = pixels[i];
        pixel_idx p = {pixels[i], idx};
        add_pixel(kdtree, p);
    }
    ptr += STARTING_SEEDS;
    int steps = 0;
    int percentage = 0;
    int reset_steps = 0;
    ptr--;
    while (++ptr < pixels + size) {
        if (++steps % (size / 100) == 0) {
            printf("Progress: %2d \tLeaves: %d\n", percentage++, kdtree->leaf_count);
            steps = 0;
        }
        if (++reset_steps % (size / 50) == 0) {
            destroy_kd_tree(kdtree);
            kdtree = build_new_tree(buffer, size_x, size_y);
            reset_steps = 0;
        }
        pixel_idx closest;
        int distance;
try_again:
        distance = find_closest(kdtree, *ptr, &closest);
        if (distance == INT_MAX) {
            printf("Somehow the tree is empty?\n");
            break;
        }
        int neighbour_idxs[8];
        int available = get_available_neighbour_idxs(closest.idx, buffer, size_x, size_y, neighbour_idxs);
        int filled_idx = -1;
        if (available) {
            i = rand_in_range(0, available - 1);
            filled_idx = neighbour_idxs[i];
            buffer[filled_idx] = *ptr;
        }
        if (available == 0) {
            remove_pixel(kdtree, closest);
            goto try_again;
        }
        if (filled_idx == -1) {
            printf("We never got to fill one in?\n");
            break;
        }
        if (available == 1) {
            /* We've exhausted this point's neighbours */
            remove_pixel(kdtree, closest);
        }
        pixel_idx new_pixel = {*ptr, filled_idx};
        /* check that we're not adding a pixel with no empty neighbours */
        available = get_available_neighbour_idxs(filled_idx, buffer, size_x, size_y, neighbour_idxs);
        if (available) {
            add_pixel(kdtree, new_pixel);
        }
    }
    destroy_kd_tree(kdtree);
}

int main() {
    /* generate all colors */
    const int no_colors =
            CHANNEL_SIZE * CHANNEL_SIZE * CHANNEL_SIZE;
    pixel* colors = malloc(no_colors * sizeof(pixel));
    pixel* ptr = colors;
    pixel* buffer = malloc(no_colors * sizeof(pixel));
    memset(buffer, 0, no_colors * sizeof(pixel));
    int r, g, b;
    int step = 256 / CHANNEL_SIZE;
    for (r = 0; r < CHANNEL_SIZE; r++) {
        for (g = 0; g < CHANNEL_SIZE; g++) {
            for (b = 0; b < CHANNEL_SIZE; b++) {
                pixel c = {r * step, g * step, b * step, 0x1};
                *ptr++ = c;
            }
        }
    }
    color_sort(colors, no_colors);
    fill(colors, buffer, IMAGE_X, IMAGE_Y);
    printf("Verifying output\n");
    for (ptr = buffer; ptr < buffer + no_colors; ptr++) {
        if (ptr->b == 0 && ptr->g == 0 && ptr->r == 0) {
            printf("Black is at %ld\n", ptr - buffer);
        }
        if (ptr->a != 0x1) {
            printf("Found unset or incorrect pixel tag value %d\n", ptr->a);
            printf("%ld\n", ptr - buffer);
            break;
        }
    }
    write_png(buffer, "test.png", IMAGE_X, IMAGE_Y);
    return 0;
}
