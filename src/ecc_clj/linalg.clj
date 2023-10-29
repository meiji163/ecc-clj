(ns ecc-clj.linalg
  (:require [ecc-clj.poly :as p]))

(def not-invertible (new Exception "matrix is not invertible"))

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

(defn swap-rows! [arr i j]
  (let [ai (aget arr i)
        aj (aget arr j)]
    (do
      (aset arr i aj)
      (aset arr j ai))))

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
      (cond
        (or (>= i rows) (>= j cols))
        mat

        (= zero (aget mat i j)) ;; find a pivot row and swap
        (let [swap-ix (first (filter
                              #(not= zero (aget mat % j))
                              (range rows)))]
          (when (nil? swap-ix) (throw not-invertible))

          (swap-rows! mat swap-ix i)
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
          ;;
          (recur (inc i) (inc j)))
        ))
    ))

(defn mat-inv [mat field]
  (let [rows (count mat)
        cols (count (first mat))]
    (if (not= rows cols)
      (throw (new Exception "matrix is not square"))

      ;; extend matrix with the identity matrix
      (let [exmat (make-array Integer/TYPE rows (* 2 cols))]
        (do
          (doall (for [i (range rows)
                       j (range cols)]
                   (aset exmat i j
                         (aget mat i j))))
          (doall (for [i (range rows)]
                   (aset exmat i (+ i rows)
                         (int 1))))

          ;; try gaussian elimination
          (let [elim (gauss-elimination exmat field)
                ;; it succeeded if the left half equals identity matrix
                success? (every? (fn [[i j]]
                                   (= (if (= i j) (int 1) (int 0))
                                      (aget elim i j)))
                                 (for [i (range rows)
                                       j (range cols)]
                                   [i j]))]

            (if (not success?) nil

                ;; extract the inverse matrix from right half
                (let [inv (make-array Integer/TYPE rows rows)]
                  (doall (for [i (range rows)
                               j (range cols)]
                           (aset inv i j
                                 (aget elim i (+ j rows)))))
                  inv))))
        ))
    ))

(defn dot [v1 v2 field]
  (let [{plus :+
         mul :*} field]
    (areduce v1 i res 0
             (plus res
                   (mul (aget v1 i)
                        (aget v2 i))))
    ))

(defn mat-vec
  "multiply matrix by vector.
  arr is an array representing a column vector."
  [mat arr field]
  (let [n (count arr)]
    (vec (for [i (range n)]
       (dot (aget mat i) arr field)))))

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
