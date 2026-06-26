/* json.h — a small, dependency-free JSON parser (first-party, for taberna).
 *
 * Parses a JSON text into a value tree and offers object/array lookup + typed
 * getters. Enough to read provenance blobs from the .mca metadata carveout
 * (e.g. roi.origin); not a full RFC 8259 validator (no \u escapes, no surrogate
 * pairs), but it handles the objects/arrays/strings/numbers/bool/null we emit.
 *
 * Ownership: json_parse() returns a heap tree freed in one call to json_free().
 * All accessor results are borrowed (owned by the tree).
 */
#ifndef TABERNA_JSON_H
#define TABERNA_JSON_H

#include <stddef.h>

typedef enum { JSON_NULL, JSON_BOOL, JSON_NUM, JSON_STR, JSON_ARR, JSON_OBJ } json_type;

typedef struct json_value json_value;
struct json_value {
  json_type type;
  double      num;     /* JSON_NUM */
  int         boolean; /* JSON_BOOL (0/1) */
  char       *str;     /* JSON_STR (NUL-terminated, decoded) */
  json_value **items;  /* JSON_ARR/JSON_OBJ element values */
  char       **keys;   /* JSON_OBJ member keys (parallel to items); NULL for arrays */
  int         count;   /* number of items */
};

/* Parse `len` bytes of JSON. Returns the root value, or NULL on syntax error. */
json_value *json_parse(const char *text, size_t len);
/* Free a tree returned by json_parse (NULL-safe). */
void json_free(json_value *v);

/* Object member by key (NULL if not an object or key absent). */
const json_value *json_obj_get(const json_value *o, const char *key);
/* Array element by index (NULL if not an array or out of range). */
const json_value *json_arr_at(const json_value *a, int i);

/* Typed reads with a default when the value is missing/wrong type. */
double      json_as_num(const json_value *v, double def);
int         json_as_int(const json_value *v, int def);
const char *json_as_str(const json_value *v, const char *def);

#endif /* TABERNA_JSON_H */
