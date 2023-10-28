(ns ecc-clj.poly-test
  (:require [ecc-clj.poly :as p]
            [ecc-clj.core :as c]
            [clojure.test :refer :all]))

(deftest poly-mod-test
  (testing "integer polynomial remainder"
    (is (= (p/mod [-4 0 -2 1] [-3 1]) [5]))
    (is (= (p/mod [-3 2 -1 1] [-2 1]) [5]))
    (is (= (p/mod [-4 4 0 5 6] [-1 1 2]) [-3 4]))
    (is (= (p/mod [0 0 0 1] [1 0 1 1]) [-1 0 -1])))
  (testing "mod2 polynomial remainder"
    (let [f p/binary-field]
      (is (= (p/mod [0 1 0 0 1] [1 1 0 1] f) [0 0 1]))
      )))

(deftest binary-mod-test
  (testing "bit-polynomial remainder"
    (is (= 4
         (p/binmod
          (c/parse-bin "11010011101100000")
          (c/parse-bin "1011")))))

  (testing "binary polynomial mod functions agree"
    (let [p1 [0 1 0 0 1]
          p2 [1 1 0 1]
          p1-bits (p/vec-to-bits p1)
          p2-bits (p/vec-to-bits p2)]
      (is (= (p/mod p1 p2 p/binary-field)
             (p/bits-to-vec
              (p/binmod p1-bits p2-bits))
             [0 0 1]))
      )))

(deftest char2-invs-test
  (testing "char2 field inverses"
    (let [prim 2
          poly (p/parse-bin "10011")
          invs (p/char2-invs prim poly)]
      (is (every? #(= 1 %)
        (map
         (fn [[a b]]
           (p/binmod (p/bin* a b) poly))
         (seq invs)))))
    ))

(deftest poly-sub-test
  (testing "integer polynomial addition"
    (is (= (p/- [0 1] [1 -1 2]) [-1 2 -2]))
    (is (= (p/+ [0 4 -2 1] [0 -4 2 2]) [0 0 0 3]))))

(deftest poly-mul-test
  (testing "integer polynomial multiplication"
    (is (= (p/* [1 -1] [0 2 1]) [0 2 -1 -1]))
    (is (= (p/* [12 34 56] [0 1]) [0 12 34 56])))
  (testing "mod2 polynomial multiplication"
    (is (= (p/* [0 1 1] [1 1 1] p/binary-field) [0 1 0 0 1]))))
