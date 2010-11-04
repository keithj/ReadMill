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

(defparameter *readmill-commands* (make-hash-table :test #'equal)
  "A hash-table mapping a ReadMill command name to a list containing
a CLI object and a function implementing the desired command.")

(defun get-cli (cmd)
  "Returns the CLI object for command string CMD."
  (first (gethash cmd *readmill-commands*)))

(defun get-fn (cmd)
  "Returns the function for command string CMD."
  (second (gethash cmd *readmill-commands*)))

(defun register-command (cmd cli fn)
  "Registers string CMD as calling FN via CLI. CLI must be a symbol
designating a CLI class."
  (check-arguments (and (stringp cmd) (symbolp cli) (functionp fn))
                   (cmd cli fn)
                   "expected a command string, a cli symbol and a function")
  (setf (gethash cmd *readmill-commands*) (list (make-instance cli) fn)))

(defun print-avail-commands (&optional (stream *error-output*))
  "Prints help for all registered commands to STREAM."
  (let ((commands (sort (loop
                           for cmd being the hash-keys of *readmill-commands*
                           collect cmd) #'string<)))
    (format stream "ReadMill version ~s~%~%" *readmill-version*)
    (write-line "Available commands:" stream)
    (terpri stream)
    (mapc (lambda (cmd)
            (let* ((cli (get-cli cmd))
                   (doc (documentation (class-name (class-of cli)) 'type)))
              (if doc
                  (help-message cli doc stream)
                  (warn "No help was found for ~a" cmd)))) commands)))

(defun readmill-cli ()
  "Applies the appropriate command line interface."
  (flet ((errmsg (c)
           (write-line (string-capitalize (format nil "~a" c) :end 1)
                       *error-output*)))
    (with-cli (argv :quit t :error-status 1)
      (let* ((cmd (first argv))
             (args (rest argv))
             (cli (get-cli cmd)))
        (handler-case
            (cond ((null cmd)
                   (print-avail-commands))
                  (cli
                   (funcall (get-fn cmd) (parse-command-line cli args) argv))
                  (t
                   (error 'unknown-command :command cmd)))
          (unknown-command (condition)
            (errmsg condition)
            (print-avail-commands)
            (quit-lisp :status 2))
          (cli-error ()
            (princ "Usage: " *error-output*)
            (let ((msg (documentation (class-name (class-of cli)) 'type)))
              (if msg
                  (help-message cli msg *error-output*)
                  (warn "No help was found for ~a~%" cmd)))
            (quit-lisp :status 3))
          (error (condition)
            (errmsg condition)
            (quit-lisp :status 4)))))))

(define-cli verbosity-mixin ()
  ((verbose "verbose" :required-option nil :value-type t
            :documentation "Print summary of action to standard output.")))

(define-cli input-file-mixin ()
  ((input-file "input-file" :required-option t :value-type 'string
               :documentation "The input file.")))

(define-cli output-file-mixin ()
  ((output-file "output-file" :required-option t :value-type 'string
                :documentation "The output file.")))

(define-cli json-log-mixin ()
  ((json-file "json-file" :required-option nil :value-type 'string
              :documentation "The JSON output file.")))

(define-cli sample-name-mixin ()
  ((sample-name "sample-name" :required-option t :value-type 'string
                :documentation "The sample name.")))

(define-cli read-group-mixin ()
  ((read-group "read-group" :required-option nil :value-type 'string
               :documentation "The read group ID.")))

(define-cli read-range-mixin ()
  ((read-start "read-start" :required-option nil :value-type 'integer
               :documentation "The zero-based start position in the read.")
   (read-end "read-end" :required-option nil :value-type 'integer
             :documentation #.(txt "The zero-based, half-open end position"
                                   "in the read."))))

(define-cli about-cli (cli)
  ((platform "platform" :required-option nil :value-type t
             :documentation "Reports the Lisp platform details.")
   (version "version" :required-option nil :value-type t
            :documentation "Reports the software version."))
  (:documentation "about [--platform] [--version]"))

(define-cli quality-plot-cli (cli read-group-mixin input-file-mixin)
  ((plot-file "plot-file" :required-option t :value-type 'string
              :documentation "The plot file."))
  (:documentation "quality-plot [--read-group <id>] --input-file <filename>
--plot-file <filename>"))

(define-cli pattern-report-cli (cli read-group-mixin input-file-mixin)
  ((report-file "report-file" :required-option t :value-type 'string
                :documentation "The plot file.")
   (pattern-char "pattern-char" :required-option t :value-type 'character
                 :documentation "The pattern character.")
   (min-freq "min-freq" :required-option t :value-type 'integer
             :documentation "The minimum pattern frequency to be reported."))
  (:documentation "pattern-report [--read-group <id>] --input-file <filename>
--report-file <filename> --pattern-char <character> --min-freq <frequency>"))

(define-cli quality-filter-cli (cli read-range-mixin
                                input-file-mixin output-file-mixin)
  ((min-quality "min-quality" :required-option t :value-type 'integer
                :documentation "The minimum quality threshold."))
  (:documentation "quality-filter --input-file <filename>
--output-file <filename> [--read-start <integer>] [--read-end <integer>]
--min-quality <integer>"))

(define-cli subseq-filter-cli (cli read-range-mixin
                               input-file-mixin output-file-mixin)
  ((queries "queries" :required-option t :value-type 'string-list
            :documentation "The subsequences to search for."))
  (:documentation "subseq-filter --input-file <filename>
--output-file <filename> [--read-start <integer>] [--read-end <integer>]
--queries <seq1,seq2 ... seqn>"))

(register-command "about" 'about-cli #'about)
(register-command "quality-plot" 'quality-plot-cli #'quality-plot)
(register-command "pattern-report" 'pattern-report-cli #'pattern-report)

(register-command "quality-filter" 'quality-filter-cli #'quality-filter)
(register-command "subseq-filter" 'subseq-filter-cli #'subseq-filter)
