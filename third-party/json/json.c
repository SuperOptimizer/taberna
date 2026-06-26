/* json.c — small recursive-descent JSON parser. See json.h. */
#include "json.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct { const char *p, *end; int ok; } jp;

static json_value *jv_new(json_type t){ json_value *v=calloc(1,sizeof *v); if(v) v->type=t; return v; }

static void skip_ws(jp *s){
  while(s->p<s->end){
    char c=*s->p;
    if(c==' '||c=='\t'||c=='\n'||c=='\r') s->p++;
    else break;
  }
}

static json_value *parse_value(jp *s);

/* Decode a JSON string starting at the opening quote; advances past the close. */
static char *parse_string_raw(jp *s){
  if(s->p>=s->end||*s->p!='"'){ s->ok=0; return NULL; }
  s->p++;
  size_t cap=16,n=0; char *out=malloc(cap);
  if(!out){ s->ok=0; return NULL; }
  while(s->p<s->end){
    char c=*s->p++;
    if(c=='"'){ out[n]='\0'; return out; }
    if(c=='\\'){
      if(s->p>=s->end) break;
      char e=*s->p++;
      switch(e){
        case '"': c='"'; break;   case '\\': c='\\'; break; case '/': c='/'; break;
        case 'n': c='\n'; break;   case 't': c='\t'; break; case 'r': c='\r'; break;
        case 'b': c='\b'; break;   case 'f': c='\f'; break;
        case 'u': /* skip 4 hex; emit '?' (we don't need unicode in provenance) */
          for(int i=0;i<4 && s->p<s->end;i++) s->p++; c='?'; break;
        default: c=e; break;
      }
    }
    if(n+1>=cap){ cap*=2; char *t=realloc(out,cap); if(!t){ free(out); s->ok=0; return NULL; } out=t; }
    out[n++]=c;
  }
  free(out); s->ok=0; return NULL;
}

static json_value *parse_object(jp *s){
  json_value *v=jv_new(JSON_OBJ); if(!v){ s->ok=0; return NULL; }
  s->p++; /* { */
  skip_ws(s);
  if(s->p<s->end && *s->p=='}'){ s->p++; return v; }
  for(;;){
    skip_ws(s);
    char *key=parse_string_raw(s);
    if(!s->ok){ free(key); json_free(v); return NULL; }
    skip_ws(s);
    if(s->p>=s->end||*s->p!=':'){ s->ok=0; free(key); json_free(v); return NULL; }
    s->p++;
    json_value *val=parse_value(s);
    if(!s->ok){ free(key); json_free(val); json_free(v); return NULL; }
    json_value **ni=realloc(v->items,(size_t)(v->count+1)*sizeof *ni);
    char **nk=realloc(v->keys,(size_t)(v->count+1)*sizeof *nk);
    if(!ni||!nk){ if(ni)v->items=ni; if(nk)v->keys=nk; free(key); json_free(val); json_free(v); s->ok=0; return NULL; }
    v->items=ni; v->keys=nk; v->items[v->count]=val; v->keys[v->count]=key; v->count++;
    skip_ws(s);
    if(s->p<s->end && *s->p==','){ s->p++; continue; }
    if(s->p<s->end && *s->p=='}'){ s->p++; return v; }
    s->ok=0; json_free(v); return NULL;
  }
}

static json_value *parse_array(jp *s){
  json_value *v=jv_new(JSON_ARR); if(!v){ s->ok=0; return NULL; }
  s->p++; /* [ */
  skip_ws(s);
  if(s->p<s->end && *s->p==']'){ s->p++; return v; }
  for(;;){
    json_value *e=parse_value(s);
    if(!s->ok){ json_free(e); json_free(v); return NULL; }
    json_value **ni=realloc(v->items,(size_t)(v->count+1)*sizeof *ni);
    if(!ni){ json_free(e); json_free(v); s->ok=0; return NULL; }
    v->items=ni; v->items[v->count++]=e;
    skip_ws(s);
    if(s->p<s->end && *s->p==','){ s->p++; skip_ws(s); continue; }
    if(s->p<s->end && *s->p==']'){ s->p++; return v; }
    s->ok=0; json_free(v); return NULL;
  }
}

static json_value *parse_value(jp *s){
  skip_ws(s);
  if(s->p>=s->end){ s->ok=0; return NULL; }
  char c=*s->p;
  if(c=='{') return parse_object(s);
  if(c=='[') return parse_array(s);
  if(c=='"'){ char *str=parse_string_raw(s); if(!s->ok)return NULL;
    json_value *v=jv_new(JSON_STR); if(!v){ free(str); s->ok=0; return NULL; } v->str=str; return v; }
  if(c=='t'||c=='f'){
    int isT=(c=='t'); const char *kw=isT?"true":"false"; size_t kl=isT?4:5;
    if((size_t)(s->end-s->p)<kl||memcmp(s->p,kw,kl)){ s->ok=0; return NULL; }
    s->p+=kl; json_value *v=jv_new(JSON_BOOL); if(!v){ s->ok=0; return NULL; } v->boolean=isT; return v;
  }
  if(c=='n'){
    if((size_t)(s->end-s->p)<4||memcmp(s->p,"null",4)){ s->ok=0; return NULL; }
    s->p+=4; return jv_new(JSON_NULL);
  }
  /* number */
  if(c=='-'||(c>='0'&&c<='9')){
    char *endp=NULL; double d=strtod(s->p,&endp);
    if(endp==s->p){ s->ok=0; return NULL; }
    s->p=endp; json_value *v=jv_new(JSON_NUM); if(!v){ s->ok=0; return NULL; } v->num=d; return v;
  }
  s->ok=0; return NULL;
}

json_value *json_parse(const char *text, size_t len){
  if(!text) return NULL;
  jp s={ .p=text, .end=text+len, .ok=1 };
  json_value *v=parse_value(&s);
  if(!s.ok){ json_free(v); return NULL; }
  return v;
}

void json_free(json_value *v){
  if(!v) return;
  if(v->type==JSON_STR) free(v->str);
  if(v->type==JSON_ARR||v->type==JSON_OBJ){
    for(int i=0;i<v->count;i++){ json_free(v->items[i]); if(v->keys) free(v->keys[i]); }
    free(v->items); free(v->keys);
  }
  free(v);
}

const json_value *json_obj_get(const json_value *o, const char *key){
  if(!o||o->type!=JSON_OBJ||!key) return NULL;
  for(int i=0;i<o->count;i++) if(o->keys[i] && !strcmp(o->keys[i],key)) return o->items[i];
  return NULL;
}
const json_value *json_arr_at(const json_value *a, int i){
  if(!a||a->type!=JSON_ARR||i<0||i>=a->count) return NULL;
  return a->items[i];
}
double json_as_num(const json_value *v, double def){ return (v&&v->type==JSON_NUM)?v->num:def; }
int    json_as_int(const json_value *v, int def){ return (v&&v->type==JSON_NUM)?(int)(v->num<0?v->num-0.5:v->num+0.5):def; }
const char *json_as_str(const json_value *v, const char *def){ return (v&&v->type==JSON_STR)?v->str:def; }
