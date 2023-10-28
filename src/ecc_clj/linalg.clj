(ns ecc-clj.linalg
  (:require [ecc-clj.poly :as p]))

(defn map-row!
  "my functional programming sin"
  [f row]
  (loop [i 0]
    (when (< i (count row))
      (aset row i
            (f (aget row i)))
      (recur (inc i)))))

(defn add-row! [r1 r2 plus]
  (loop [i 0]
    (when (< i (count r1))
      (aset r1 i
            (plus (aget r1 i) (aget r2 i)))
      (recur (inc i)))))

(defn gauss-elimination
  [mat field]
  (let [{plus :+
         mul :*
         div :/
         sub :-
         zero :zero} field

        rows (count mat)
        cols (count (first mat))]
    (loop [i 0
           j 0]
      (let [_ (clojure.pprint/pprint mat)
            _ (println "i:" i " j:" j)]
        (cond
          (or (>= i rows) (>= j cols))
          mat

          (= zero (aget mat i j)) ;; find a pivot row and swap
          (let [swap-ix (first (filter
                                #(not= zero (aget % j))
                                (range cols)))]
            (swap-rows! arr swap-ix i)
            (recur i j))

          :else
          ;; subtract current row from other rows
          ;; such that current column has all zeros
          (do
            (doall
             (map
              (fn [k] (when (not= i k)
                        (let [aij (aget mat i j)
                              akj (aget mat k j)
                              ai (aget mat i)
                              factor (if (= zero akj)
                                       0
                                       (div akj aij))
                              sub-row (amap ai i ret
                                            (mul factor (aget ai i)))]
                          (add-row! (aget mat k) sub-row plus))))
              (range rows)))
            ;; make this diagonal = 1
            (let [aij (aget mat i j)
                  ai (aget mat i)]
              (map-row! #(div % aij) ai))

            (recur (inc i) (inc j)))
          )))
    ))

(defn matrix* [m1 m2 field]
  (let [{plus :+
         mul :*} field
        rows1 (count m1)
        cols1 (count (first m1))
        rows2 (count m2)
        cols2 (count (first m2))
        rows3 (count m1)
        res (make-array Integer/TYPE rows1 cols2)]
    (cond
      (or (not= cols1 rows2)) (throw (new Exception "incorrect matrix sizes"))
      :else
      "TODO")
    ))

(def GF8
  (let [GF2-poly (Integer/parseInt "1011" 2)]
    (p/char2-field 2 GF2-poly)))

(println "TESTTTT")
(def arr (make-array Integer/TYPE 2 4))
(aset arr 0 0 (int 7))
(aset arr 0 1 (int 3))
(aset arr 1 1 (int 2))
(aset arr 0 2 (int 1))
(aset arr 1 3 (int 1))


(clojure.pprint/pprint arr)

(gauss-elimination arr GF8)

(clojure.pprint/pprint arr)
