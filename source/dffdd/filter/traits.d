module dffdd.filter.traits;


enum bool isUpdater(T, size_t P, C) = is(typeof((T t){
    C[P][] w, x;
    C x00, e;
    t.update(w, x, x00, e);
}));
