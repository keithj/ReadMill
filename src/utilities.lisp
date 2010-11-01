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

(defvar *software-name* "ReadMill")
(defvar *software-version* "0.0.3"
  "The version number of the ReadMill application.")

(deftype quality-score ()
  '(integer 0 100))

(defun peek-alignment (bam)
  "Returns the next alignment in BAM stream, or NIL."
  (let ((pos (bgzf-tell bam)))
    (prog1
        (read-alignment bam)
      (bgzf-seek bam pos))))

(defun peek-read-length (bam)
  "Returns the length of the next read BAM stream."
  (let ((aln (peek-alignment bam)))
    (if aln
        (read-length aln)
        (error 'invalid-operation-error
               :format-control (txt "cannot peek read length;"
                                    "reached end of stream")))))

(defun decode-phred-quality (c)
  "Returns the Phred quality score encoded by the character C."
  (- (char-code c) 33))

(defun vector-mean (sums count)
  (flet ((mean (sum)
           (float (/ sum count))))
    (let ((means (make-array (length sums) :element-type 'float
                             :initial-element 0.f0)))
      (map-into means #'mean sums))))

(defun plot-per-base (plot-filespec values sample-name y-label)
  (plot-png nil values plot-filespec :title sample-name
            :x-label "Base position" :y-label y-label))

(defun plot-png (x y png-filespec &key title x-label y-label)
  "Plots X against Y and saves the result to file PNG-FILESPEC. The
plot is annotated with TITLE, X-LABEL and Y-LABEL."
  (let ((plotter (dxr:run-gnuplot))
        (plot (make-instance
               'dxr:2d-plot
               :title title
               :x-axis (make-instance 'dxr:axis :label x-label :position :x)
               :y-axis (make-instance 'dxr:axis :label y-label :position :y)
               :series (make-instance 'dxr:xy-series
                                      :x-values (or x (iota (length y)))
                                      :y-values y
                                      :style '(:linespoints
                                               :smooth
                                               :csplines)))))
    (dxr:draw-plot plotter plot :terminal :png :output png-filespec)
    (dxr:stop-gnuplot plotter)))
