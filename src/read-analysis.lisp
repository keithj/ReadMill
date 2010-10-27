;;;
;;; Copyright (C) 2010 Genome Research Ltd. All rights reserved.
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

(in-package :uk.ac.sanger.readmill)

(defun plot-mean-qualities (plot-filespec sample-name bam-filespec)
  "Writes a mean base quality plot for BAM file BAM-FILESPEC to
PLOT-FILESPEC using the plot title SAMPLE-NAME. Will accept reads of
mixed length without warning."
  (with-bgzf (bam bam-filespec)
    (do* ((aln (read-alignment bam) (read-alignment bam))
          (rlen 0 (if aln
                      (read-length aln)
                      rlen))
          (sums (make-array rlen :element-type 'fixnum :initial-element 0)
                (if (> rlen (length sums))
                    (replace (make-array rlen :element-type 'fixnum
                                         :initial-element 0) sums)
                    sums))
          (count 0 (1+ count)))
         ((null aln) (plot-per-base plot-filespec (vector-mean sums count)
                                    sample-name "Mean base quality"))
      (loop
         with qual = (quality-string aln)
         for i from 0 below rlen
         do (incf (aref sums i) (decode-phred-quality (char qual i)))))))

(defun base-patterns (char mask-char bam-filespec &optional (min-freq 1))
  "Returns a list of strings and integers describing the most common
patterns of CHAR in the read sequences from BAM-FILESPEC and their
frequencies. CHAR must occur at least MIN-FREQ times in the sequence
to be reported. All other characters in the sequences are masked with
MASK-CHAR."
  (flet ((make-mask (len positions)
           (let ((mask (make-array len :element-type 'base-char
                                  :initial-element mask-char)))
             (loop
                for i in positions
                do (setf (char mask i) char)
                finally (return mask)))))
    (with-bgzf (bam bam-filespec)
      (read-bam-meta bam)
      (loop
         with pos-table = (make-hash-table :test #'equal)
         for aln = (read-alignment bam)
         while aln
         count aln into read-count
         maximize (read-length aln) into read-length
         do (let ((pos (vector-positions char (seq-string aln) :test #'char=)))
              (when (>= (length pos) min-freq)
                (if (gethash pos pos-table)
                    (incf (gethash pos pos-table))
                    (setf (gethash pos pos-table) 1))))
         finally (return
                   (loop
                      for pos being the hash-keys in pos-table
                      using (hash-value freq)
                      collect (list freq (make-mask read-length pos))
                      into summary
                      finally (return (values summary read-count))))))))
