(ns ecc-clj.core-test
  (:require [clojure.test :refer :all]
            [ecc-clj.core :refer :all]))

(deftest mod-inverses
  (testing "mod inverses"
    (let [p 17]
      (is
       (every? #(= 1 %)
               (map
                (fn [[a b]] (mod (* a b) p))
                (seq (mod-invs p))))))))

(deftest poly-mod-test
  (testing "integer polynomial remainder"
    (is (= (poly-mod [-4 0 -2 1] [-3 1]) [5]))
    (is (= (poly-mod [-3 2 -1 1] [-2 1]) [5]))
    (is (= (poly-mod [-4 4 0 5 6] [-1 1 2]) [-3 4]))
    (is (= (poly-mod [0 0 0 1] [1 0 1 1]) [-1 0 -1])))
  (testing "mod2 polynomial remainder"
    (let [f binary-field]
      (is (= (poly-mod [0 1 0 0 1] [1 1 0 1] f) [0 0 1]))
      )))

(deftest bpoly-mod-test
  (testing "bit-polynomial remainder"
    (is (= 4
         (bpoly-mod
          (Integer/parseInt "11010011101100000" 2)
          (Integer/parseInt "1011" 2)))))

  (testing "binary polynomial mod functions agree"
    (let [p1 [0 1 0 0 1]
          p2 [1 1 0 1]
          p1-bits (vec-to-bits p1)
          p2-bits (vec-to-bits p2)]
      (is (= (poly-mod p1 p2 binary-field)
             (bits-to-vec
              (bpoly-mod p1-bits p2-bits))
             [0 0 1]))
      )))

(deftest char2-invs-test
  (testing "char2 field inverses"
    (let [prim 2
          poly (parse-bin "10011")
          invs (char2-invs prim poly)]
      (is (every? #(= 1 %)
        (map
         (fn [[a b]]
           (bpoly-mod (bpoly-mul a b) poly))
         (seq invs)))))
    ))

(deftest poly-sub-test
  (testing "integer polynomial addition"
    (is (= (poly-sub [0 1] [1 -1 2]) [-1 2 -2]))
    (is (= (poly-add [0 4 -2 1] [0 -4 2 2]) [0 0 0 3]))))

(deftest poly-mul-test
  (testing "integer polynomial multiplication"
    (is (= (poly-mul [1 -1] [0 2 1]) [0 2 -1 -1]))
    (is (= (poly-mul [12 34 56] [0 1]) [0 12 34 56])))
  (testing "mod2 polynomial multiplication"
    (is (= (poly-mul [0 1 1] [1 1 1] binary-field) [0 1 0 0 1]))))
