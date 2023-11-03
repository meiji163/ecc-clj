(ns ecc-clj.core-test
  (:require [clojure.test :refer :all]
            [ecc-clj.core :refer :all]
            [ecc-clj.poly :as p]))

(deftest mod-inverses
  (testing "mod inverses"
    (let [p 17]
      (is
       (every? #(= 1 %)
               (map
                (fn [[a b]] (mod (* a b) p))
                (seq (p/mod-invs p))))))))

(deftest RS-7-5-test
  (testing "encoding"
    (let [data [5 2 3 1 6]]
      (is (= (encode data RS-7-5 GF8)
             [6 6 5 2 3 1 6]))))
  (testing "decoding without errors"
    (let [data [0 4 3 0 1]]
      (is (= data
             (-> data
                 (encode RS-7-5 GF8)
                 (RS-7-5-decode))))))
  (testing "decoding with errors"
    (let [encoded-poly [6 6 5 2 3 1 6]
          data-poly [5 2 3 1 6]]
      (is
       (every? #(= data-poly %)
        (for [err (single-errors 7 8)]
          (let [err-poly (p/+ encoded-poly err GF8)]
            (RS-7-5-decode err-poly)))))
      )))

(deftest RS-7-3-test
  (testing "encoding"
    (let [data [1 6 3]]
      (is (= (encode data RS-7-3 GF8)
             [3 3 2 6 1 6 3])))))
