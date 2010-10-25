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

(define-condition readmill-error (error simple-text-condition)
  ()
  (:report (lambda (condition stream)
             (format stream "ReadMill error~@[: ~a~]"
                     (message-of condition))))
  (:documentation "An error that is raised when using ReadMill."))

(define-condition readmill-warning (warning simple-text-condition)
  ()
  (:report (lambda (condition stream)
             (format stream "ReadMill warning~@[: ~a~]"
                     (message-of condition))))
  (:documentation "An warning that is raised when using ReadMill."))
