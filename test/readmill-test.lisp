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

(defmacro with-fake-bam ((filespec basename &rest args) &body body)
  "Evaluates BODY with FILESPEC bound to the file name of a newly
created fake BAM file defined by ARGS. If no error occurs in BODY, the
fake BAM file will be deleted automatically."
  `(let ((,filespec (fake-bam-file (dxi:tmp-pathname
                                    :basename ,basename
                                    :tmpdir (merge-pathnames "data")
                                    :type "bam")
                                   ,@args)))
     (prog1
         (progn
           ,@body)
       (delete-file ,filespec))))       ; do not delete on error in body

(addtest (readmill-tests) base-patterns/1
  (with-fake-bam (bamfile "base-patterns"
                          :num-refs 1 :ref-length 500 :read-length 10
                          :seq-fn (lambda (length)
                                    (ambiguous-positions length 5 7 9)))
    (with-bam (in () bamfile)
      (let ((expected '((".....N.N.N" . 10)))
            (patterns (base-patterns in #\N #\.)))
        (ensure (equalp expected patterns)
                :report "Expected ~a but found ~a"
                :arguments (expected patterns))))))

(addtest (readmill-tests) quality-filter/1
  (with-fake-bam (bamfile "quality-filter"
                          :num-refs 1 :ref-length 500 :read-length 10
                          :qual-fn (lambda (length)
                                     (quality-ranges length 20 10 0 6)))
    (let ((filter1 (readmill::make-quality-p 20 :start 0 :end 6))
          (filter2 (make-quality-p 20 :start 6 :end 10)))
      (with-bam (in () bamfile)         ; filter 10/10
        (let* ((fin (discarding-if filter1 in))
               (n (loop
                     while (has-more-p fin)
                     count (next fin))))
          (ensure (zerop n)
                  :report "Expected 0 but found ~d"
                  :arguments (n))))
      (with-bam (in () bamfile)         ; filter 0/10
        (let* ((fin (discarding-if filter2 in))
               (n (loop
                     while (has-more-p fin)
                     count (next fin))))
          (ensure (= 10 n)
                  :report "Expected 10 but found ~d"
                  :arguments (n)))))))

(addtest (readmill-tests) read-group-filter/1
  (flet ((read-group (aln)
           (assocdr :rg (alignment-tag-values aln))))
    (with-fake-bam (bamfile "read-group-filter"
                            :num-refs 2 :ref-length 500 :read-length 10)
      (let ((filter (make-rg-p "0")))
        (with-bam (in () bamfile)
          (let* ((fin (discarding-if filter in))
                 (n (loop
                       while (has-more-p fin)
                       count (string= "0" (read-group (next fin))))))
            (ensure (= 10 n)
                    :report "Expected 10 but found ~d"
                    :arguments (n))))))))

(addtest (readmill-tests) subseq-filter/1
  (with-fake-bam (bamfile "subseq-filter"
                          :num-refs 1 :ref-length 500 :read-length 10
                          :seq-fn (lambda (length)
                                    (common-subsequence length "NNNN")))
    (let ((filter1 (make-subseq-p "NNNN" :start 0 :end 4))
          (filter2 (make-subseq-p "NNNN" :start 4 :end 10)))
      (with-bam (in () bamfile)
        (let* ((fin (discarding-if filter1 in))
               (n (loop
                     while (has-more-p fin)
                     count (next fin))))
          (ensure (zerop n)
                  :report "Expected 0 but found ~d"
                  :arguments (n))))
      (with-bam (in () bamfile)
        (let* ((fin (discarding-if filter2 in))
               (n (loop
                     while (has-more-p fin)
                     count (next fin))))
          (ensure (= 10 n)
                  :report "Expected 10 but found ~d"
                  :arguments (n)))))))
