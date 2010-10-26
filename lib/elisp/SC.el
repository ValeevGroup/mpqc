;;
;; SC.el: stuff for formatting files in the SC libraries
;;
;; Copyright (C) 1996 Limit Point Systems, Inc.
;;
;; Author: Curtis Janssen <cljanss@ca.sandia.gov>
;; Maintainer: SNL
;;
;; This file is part of MPQC.
;;
;; MPQC is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.
;;
;; MPQC is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;;
;; You should have received a copy of the GNU General Public License
;; along with the MPQC; see the file COPYING.  If not, write to
;; the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
;;
;; The U.S. Government is granted a limited license as per AL 91-7.
;;

(require 'cc-mode)
(cond ((> emacs-major-version 19) (c-initialize-cc-mode)))

(setq clj-c-basic-half-offset 2)

(defun clj-adaptive-block-open (langelem)
  ;; when substatement is on semantics list, return
  ;; -(c-basic-offset - clj-c-basic-half-offset) to give a
  ;; total offset of clj-c-basic-half-offset,
  ;; otherwise return clj-c-basic-half-offset
  (if (assq 'substatement c-semantics)
      (+ clj-c-basic-half-offset (- c-basic-offset))
    clj-c-basic-half-offset))

(defun clj-lineup-math (langelem)
  ;; line up math statement-cont so that stuff after the "+", "-", etc
  ;; lines up with the stuff after the equals
  (save-excursion
    (let ((adjustment (progn
                        (beginning-of-line)
                        (skip-chars-forward " \t" (c-point 'eol))
                        (- (current-column)
                           (progn (skip-chars-forward " \t+-/*" (c-point 'eol))
                                  (current-column)))))
          (curcol (progn
		    (goto-char (cdr langelem))
		    (current-column))))
      (skip-chars-forward "^=" (c-point 'eol))
      (if (/= (following-char) ?=)
	  ;; there's no equal sign on the line
	  c-basic-offset
	;; calculate indentation column after equals and ws and sign
	(forward-char 1)
	(skip-chars-forward " \t-")
	(+ (- (current-column) curcol) adjustment))
      )))

(defun clj-adaptive-block-close (langelem)
  ;; these closes blocks in a way that is consistent with the way
  ;; clj-adaptive-statement-block-intro indents the first statement
  (- (clj-adaptive-statement-block-intro langelem)
     (- c-basic-offset clj-c-basic-half-offset))
)

(defun clj-adaptive-statement-block-intro (langelem)
  ;; this lines up the first statement in a block by a full basic
  ;; offset, unless we are lining up to a "{" which is already
  ;; half indented
  (save-excursion
    (progn
      (goto-char (cdr langelem))
      (if (/= (following-char) ?{)
          ;; next char is not a "{"
          c-basic-offset
        ;; use remainder of half offset
        (- c-basic-offset clj-c-basic-half-offset))
      )))

(defun clj-condensed-adaptive-statement-block-intro (langelem)
  ;; this lines up the first statement in a block by a full basic
  ;; offset, unless we are lining up to a "{" which is already
  ;; indented
  (save-excursion
    (progn
      (goto-char (cdr langelem))
      (if (/= (following-char) ?{)
          ;; next char is not a "{"
          c-basic-offset
        ;; we're already indendted
        0)
      )))

(defun clj-condensed-adaptive-block-close (langelem)
  ;; these closes blocks in a way that is consistent with the way
  ;; clj-condensed-adaptive-statement-block-intro indents the first statement
  (clj-condensed-adaptive-statement-block-intro langelem)
)

;;
;; this is the style to use when editting Ed's files
;;
(c-add-style "ETS" '((c-basic-offset . 2)
                     (c-offsets-alist . ((access-label      . -)
                                         (inclass           . ++)
                                         (label             . 0)
                                         ))
                     )
             )

;;
;; this is the style to use when editing Curt's files
;;
(c-add-style "CLJ" '(
    (c-offsets-alist . (
        (block-open      . clj-adaptive-block-open)
        (statement       . c-lineup-runin-statements)
        (statement-cont  . clj-lineup-math)
        (statement-block-intro . clj-adaptive-statement-block-intro)
        (defun-block-intro . 2)
        (inher-intro . 2)
        (access-label . -2)
        (block-close . clj-adaptive-block-close)
        (member-init-intro . 2)
        )
    ))
)

;;
;; Curt's other style
;;
(c-add-style "CLJ-CONDENSED" '(
    ;(c-echo-syntactic-information-p . t)
    (c-basic-offset . 2)
    (c-offsets-alist . (
       (statement-block-intro . clj-condensed-adaptive-statement-block-intro)
       (statement-cont . c-lineup-math)
       (inclass . ++)
       (access-label . -)
       (block-close . clj-condensed-adaptive-statement-block-intro)
       (substatement-open . +)
       (block-open . +)
       )
    ))
)

(defun clj-condensed-style ()
  "Change to condensed C indentation"
  (interactive)
  (c-set-style "CLJ-CONDENSED")
  )
(defun clj-style ()
  "Change to insane C indentation"
  (interactive)
  (c-set-style "CLJ")
  )
(defun ets-style ()
  "Change to sensible C indentation"
  (interactive)
  (c-set-style "ETS")
  )

(define-key c-mode-map "\C-ce" 'ets-style)
(define-key c-mode-map "\C-cj" 'clj-style)
(define-key c-mode-map "\C-cc" 'clj-condensed-style)
(define-key c-mode-map "\C-j"  'reindent-then-newline-and-indent)
(define-key c-mode-map "\C-m"  'newline-and-indent)

(define-key c++-mode-map "\C-ce" 'ets-style)
(define-key c++-mode-map "\C-cj" 'clj-style)
(define-key c++-mode-map "\C-cc" 'clj-condensed-style)
(define-key c++-mode-map "\C-j"  'reindent-then-newline-and-indent)
(define-key c++-mode-map "\C-m"  'newline-and-indent)

(define-key java-mode-map "\C-ce" 'ets-style)
(define-key java-mode-map "\C-cj" 'clj-style)
(define-key java-mode-map "\C-cc" 'clj-condensed-style)
(define-key java-mode-map "\C-j"  'reindent-then-newline-and-indent)
(define-key java-mode-map "\C-m"  'newline-and-indent)

;;
;; stuff for CLJ's compile hacks
;;

(defun compile-modify-path (thisdir)
  (let ((tmpdir (expand-file-name thisdir)))
    (setq thisdir "")
    (while (>= (length tmpdir) (length sc-src-dir))
      (if (string= (substring tmpdir 0 (length sc-src-dir)) sc-src-dir)
          (let ()
            (setq thisdir (concat thisdir sc-arch-dir))
            (setq tmpdir (substring tmpdir (length sc-src-dir) nil))
            )
        (let ()
          (setq thisdir (concat thisdir (substring tmpdir 0 1)))
          (setq tmpdir (substring tmpdir 1 nil))
          )
        )
      )
    (setq thisdir (concat thisdir tmpdir))
    )
  thisdir
)

;;
;; stuff for inserting copyleft notices
;;

(defvar copyleft-owner "Limit Point Systems, Inc."
  "This is the owner of the copyleft.  Defaults to LPS.")

(defun set-copyleft-owner (owner)
  "Set the copyleft-owner variable."
  (interactive (list (read-from-minibuffer "Copyleft Owner: "
                                           copyleft-owner nil nil nil)))
  (setq copyleft-owner owner))

(defvar copyleft-author user-full-name
  "This is the author of the file.  Defaults to the user editing the file.")

(defun set-copyleft-author (author)
  "Set the copyleft-author variable."
  (interactive (list (read-from-minibuffer "Author: "
                                           copyleft-author nil nil nil)))
  (setq copyleft-author author))

(defvar copyleft-address user-mail-address
  "This is the email address of the author of the file.  Defaults to the
   address of the user editing the file.")

(defun set-copyleft-address (address)
  "Set the copyleft-address variable."
  (interactive (list (read-from-minibuffer "E-mail address: "
                                           copyleft-address nil nil nil)))
  (setq copyleft-address address))

(defvar copyleft-maintainer "LPS"
  "This is the official maintaner of the file. Defaults to LPS")

(defun set-copyleft-maintainer (maintainer)
  "Set the copyleft-maintainer variable."
  (interactive (list (read-from-minibuffer "Maintainer: "
                                           copyleft-maintainer nil nil nil)))
  (setq copyleft-maintainer maintainer))

(defvar copyleft-default-comment-start "#"
  "The default symbol to use to begin a comment. Defaults to \"#\".")

(defvar copyleft-default-comment-cont "#"
  "The default symbol to use to continue a comment. Defaults to \"#\".")

(defvar copyleft-default-comment-end "#"
  "The default symbol to use to end a comment. Defaults to \"#\".")

(defun copyleft-set-comments (start cont end)
  "Set the comment symbols.

   (copyleft-set-comments START CONT END)"
  (interactive (list (read-from-minibuffer "Comment start: "
                              copyleft-default-comment-start nil nil nil)
                     (read-from-minibuffer "Comment continue: "
                              copyleft-default-comment-cont nil nil nil)
                     (read-from-minibuffer "Comment end: "
                              copyleft-default-comment-end nil nil nil)
                     ))

  (setq copyleft-default-comment-start start)
  (setq copyleft-default-comment-cont cont)
  (setq copyleft-default-comment-end end)
)

(defun insert-copyleft ()
  "Insert the notice."
  (interactive)
  (set-window-point (display-buffer (current-buffer)) (point-min))
  (cond ((eq major-mode 'c++-mode)
         (setq comment-start "//")
         (setq comment-cont "//")
         (setq comment-end "//")
         )
        ((eq major-mode 'c-mode)
         (setq comment-start "/*")
         (setq comment-cont " *")
         (setq comment-end " */")
         )
        ((eq major-mode 'emacs-lisp-mode)
         (setq comment-start ";;")
         (setq comment-cont ";;")
         (setq comment-end ";;")
         )
        ((eq major-mode 'makefile-mode)
         (setq comment-start "#")
         (setq comment-cont "#")
         (setq comment-end "#")
         )
        ('t
         (setq comment-start copyleft-default-comment-start)
         (setq comment-cont copyleft-default-comment-cont)
         (setq comment-end copyleft-default-comment-end)
         )
        )

  (setq description (read-from-minibuffer "Description: "))

  (insert comment-start "\n")
  (insert comment-cont " " (file-name-nondirectory buffer-file-name) 
      (cond ((not (string= description "")) (concat " --- " description "\n"))
            ('t "\n")))
  (insert comment-cont "\n")
  (insert comment-cont " Copyright (C) "
          (substring (current-time-string) 20) " "
          copyleft-owner "\n")
  (insert comment-cont "\n")
  (insert comment-cont " Author: " copyleft-author " <" copyleft-address ">\n")
  (insert comment-cont " Maintainer: " copyleft-maintainer "\n")
  (insert comment-cont "\n")
  (insert comment-cont " This file is part of the SC Toolkit.\n")
  (insert comment-cont "\n")

  (insert comment-cont " The SC Toolkit is free software; you can redistribute it and/or modify\n")
  (insert comment-cont " it under the terms of the GNU Library General Public License as published by\n")
  (insert comment-cont " the Free Software Foundation; either version 2, or (at your option)\n")
  (insert comment-cont " any later version.\n")
  (insert comment-cont "\n")

  (insert comment-cont " The SC Toolkit is distributed in the hope that it will be useful,\n")
  (insert comment-cont " but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
  (insert comment-cont " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n")
  (insert comment-cont " GNU Library General Public License for more details.\n")
  (insert comment-cont "\n")

  (insert comment-cont " You should have received a copy of the GNU Library General Public License\n")
  (insert comment-cont " along with the SC Toolkit; see the file COPYING.LIB.  If not, write to\n")
  (insert comment-cont " the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.\n")
  (insert comment-cont "\n")

  (insert comment-cont " The U.S. Government is granted a limited license as per AL 91-7.\n")
  (insert comment-end "\n")
)

(define-key c-mode-map "\C-ci" 'insert-copyleft)
(define-key c++-mode-map "\C-ci" 'insert-copyleft)
(define-key java-mode-map "\C-ci" 'insert-copyleft)
