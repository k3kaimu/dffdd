module dffdd.filter.traits;


enum bool isUpdater(T) = is(typeof((T t){
    float[] w, x;
    float e;
    t.start(w);
    t.update(w, x, e);
}));
