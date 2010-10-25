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
    (read-bam-meta bam)
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
