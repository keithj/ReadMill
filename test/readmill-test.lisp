;;
;;; Copyright (C) 2008-2010 Genome Research Ltd. All rights reserved.
;;;
;;; This file is part of readmill.
;;;
;;; This program is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;;

(in-package :uk.ac.sanger.readmill-test)

(deftestsuite readmill-tests ()
  ())

(make-random-state)
(defun random-base ()
  (ecase (random 4)
    (0 #\A)
    (1 #\C)
    (2 #\G)
    (3 #\T)))

(defun ambiguous-positions (length &rest positions)
  "Returns a string of LENGTH characters with Ns at POSITIONS and
random bases elsewhere."
  (loop
     with seq = (make-string length)
     for i from 0 below length
     for base = (if (member i positions)
                    #\N
                    (random-base))
     do (setf (char seq i) base)
     finally (return seq)))

(defun quality-ranges (length quality1 quality2 &optional (start 0) end)
  "Returns a string of LENGTH characters with encoded Phred scores of
QUALITY1 or QUALITY2 between START and END."
  (let ((end (or end length)))
    (loop
       with qual = (make-string length :initial-element
                                (encode-phred-quality quality1))
       for i from start below end
       do (setf (char qual i) (encode-phred-quality quality2))
       finally (return qual))))

(defun common-subsequence (length subseq &optional (position 0))
  "Returns a string of LENGTH characters with string SUBSEQ at
POSITION and random bases elsewhere."
  (loop
     with seq = (make-string length)
     for i from 0 below length
     do (setf (char seq i) (random-base))
     finally (return (replace seq subseq :start1 position))))

(defun encode-phred-quality (q)
  "Returns the character encoding Phred quality Q."
  (code-char (+ q 33)))

(defmacro with-tmp-test-file ((file &optional (basename "") type) &body body)
  "Sets up a temporary pathname in ./data for testing. Deletes any
file created, unless a test-condition is raised."
  `(handler-bind ((test-condition #'delete-tmp-pathname))
     (with-tmp-pathname (,file :tmpdir (merge-pathnames "data")
                               :basename ,basename :type ,type)
       ,@body)))

(defun test-count (file expected &optional filter)
  (with-bam (in () file)
    (let* ((fin (if filter
                    (discarding-if filter in)
                    in))
           (n (loop
                 while (has-more-p fin)
                 count (next fin))))
      (ensure (= expected n)
              :report "Expected ~d but found ~d"
              :arguments (expected n)))))

(addtest (readmill-tests) base-patterns/1
  (with-tmp-test-file (bam-file "base-patterns-" "bam")
    (let* ((rg "read_group_0")
           (num-refs 1)
           (ref-len 500)
           (gen (alignment-generator
                 0 rg :read-length 10 :insert-length 380 :step-size 10
                 :end ref-len :seq-fn (lambda (length)
                                        (ambiguous-positions length 5 7 9))))
           (bam (generate-bam-file bam-file num-refs ref-len (list rg) gen)))
      (with-bam (in () bam)
        (let ((expected '((".....N.N.N" . 20)))
              (patterns (base-patterns in #\N #\.)))
          (ensure (equalp expected patterns)
                  :report "Expected ~a but found ~a"
                  :arguments (expected patterns)))))))

(addtest (readmill-tests) quality-filter/1
  (with-tmp-test-file (bam-file "quality-filter-" "bam")
    (let* ((rg "read_group_0")
           (num-refs 1)
           (ref-len 500)
           (gen (alignment-generator
                 0 rg :read-length 10 :insert-length 380 :step-size 10
                 :end ref-len :quality-fn (lambda (length)
                                            (quality-ranges length 20 10 0 6))))
           (bam (generate-bam-file bam-file num-refs ref-len (list rg) gen))
           (filter1 (make-quality-p 20 :start 0 :end 6))
           (filter2 (make-quality-p 20 :start 6 :end 10)))
      (test-count bam 0 filter1)       ; filter 20/20
      (test-count bam 20 filter2))))   ; filter 0/20

(addtest (readmill-tests) read-group-filter/1
  (flet ((read-group (aln)
           (assocdr :rg (alignment-tag-values aln))))
    (with-tmp-test-file (bam-file "read-group-filter-" "bam")
      (let* ((rg0 "read_group_0")
             (rg1 "read_group_1")
             (num-refs 2)
             (ref-len 500)
             (gen1 (alignment-generator
                    0 rg0 :read-length 10 :insert-length 380 :step-size 10
                    :end ref-len))
             (gen2 (alignment-generator
                    0 rg1 :read-length 10 :insert-length 380 :step-size 10
                    :end ref-len))
             (bam (generate-bam-file bam-file num-refs ref-len (list rg0 rg1)
                                     gen1 gen2)))
        (with-bam (in () bam)
          (let* ((fin (discarding-if (complement (make-rg-p rg0)) in))
                 (expected 20)
                 (n (loop
                       while (has-more-p fin)
                       count (string= rg0 (read-group (next fin))))))
            (ensure (= expected n)
                    :report "Expected ~d but found ~d"
                    :arguments (expected n))))))))

(addtest (readmill-tests) subseq-filter/1
  (with-tmp-test-file (bam-file "subseq-filter-" "bam")
    (let* ((rg "read_group_0")
           (num-refs 1)
           (ref-len 500)
           (gen (alignment-generator
                 0 rg :read-length 10 :insert-length 380 :step-size 10
                 :end ref-len :seq-fn (lambda (length)
                                        (common-subsequence length "NNNN"))))
           (bam (generate-bam-file bam-file num-refs ref-len (list rg) gen))
           (filter1 (make-subseq-p "NNNN" :start 0 :end 4))
           (filter2 (make-subseq-p "NNNN" :start 4 :end 10)))
      (test-count bam 0 filter1 )       ; filter 20/20
      (test-count bam 20 filter2))))    ; filter 0/20

(addtest (readmill-tests) orphan-filter/1
  (with-tmp-test-file (sorted-file "name-sorted-")
    (with-tmp-test-file (in-file "orphan-filter-in-")
      ;; Make a BAM file containing orphans
      (let* ((rg "read_group_0")
             (num-refs 1)
             (ref-len 500)
             (gen (alignment-orphanizer
                   (alignment-generator
                    0 rg :read-length 10 :insert-length 380 :step-size 10
                    :end ref-len)))
             (in-file (generate-bam-file
                       in-file num-refs ref-len (list rg) gen)))
        (test-count in-file 15))        ; dropped alternate last frags
      ;; Must be name sorted to detect orphans
      (sort-bam-file in-file sorted-file :sort-order :queryname))
    ;; Rewrite the BAM file without orphans
    (with-tmp-test-file (out-file "orphan-filter-out-" "bam")
      (with-bam (in (header num-refs ref-meta) sorted-file)
        (with-bam (out (header num-refs ref-meta) out-file
                       :direction :output :if-does-not-exist :create)
          (let ((out (pair-consumer out)))
            (loop
               while (has-more-p in)
               do (consume out (next in))))))
      (test-count out-file 10)))) ; dropped alternate last frags and mates

(addtest (readmill-tests) split-bam/1
  (with-tmp-test-file (in-file "split-in-" "bam")
    (let* ((rg "read_group_0")
           (num-refs 1)
           (ref-len 500)
           (gen (alignment-generator
                 0 rg :read-length 10 :insert-length 380 :step-size 10
                 :end ref-len))
           (in-file (generate-bam-file
                     in-file num-refs ref-len (list rg) gen)))
      ;; File contains 20 reads
      (with-tmp-test-file (out-name "split-out-")
        (multiple-value-bind (parts sizes)
            (split-bam-by-size (list "test") in-file out-name :max-size 9)
          (let ((expected (mapcar (lambda (n)
                                    (format nil "~a.part.~d.bam"
                                            (file-namestring out-name) n))
                                  '(0 1 2)))
                (found (mapcar #'file-namestring parts)))
          (ensure (equalp expected found)
                  :report "Expected ~a but found ~a"
                  :arguments (expected found))
          (ensure (equal '(9 9 2) sizes))
          (mapcar #'delete-file parts)))))))
