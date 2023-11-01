(ns ecc-clj.linalg-test
  (:require [ecc-clj.linalg :as la]
            [ecc-clj.core :as c]
            [clojure.test :refer :all]))

(deftest mat-inv-test
  (testing "invert 2x2 GF8 matrix"
    (let [mat (into-array
               [(into-array Integer/TYPE [7 3])
                (into-array Integer/TYPE [0 2])])
          inv (la/mat-inv mat c/GF8)]
      (is (= [4 6] (vec (aget inv 0))))
      (is (= [0 5] (vec (aget inv 1)))))))
