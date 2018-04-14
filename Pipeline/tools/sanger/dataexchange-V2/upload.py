#!/usr/bin/env python
#-*- coding: utf-8 -*-
import argparse
import os
from upload_task import UploadTask
from get_file_list import get_list


parser = argparse.ArgumentParser(description="根据验证码，上传一个文件夹至某一个项目下")
parser.add_argument("-i", "--input_dir", help="输入的文件夹的路径", required=True)
parser.add_argument("-l", "--file_list", help="文件列表名字,不存在则根据input_dir自动生成, 用于描述文件的信息, 其表头必须为文件路径、文件大小、文件描述、是否锁定，锁定文件不允许客户类型用户下载", required=True)
# parser.add_argument("-g", "--get_list", help="默认F，如果为T：则根据input_dir生成一个新文件列表文件,如果为F：则自动使用list文件或根据inputdir生成list文件后直接开始上传.",default="F")

parser.add_argument("-c", "--identity_code", help="验证码,上传时必须提供", default="none")
parser.add_argument("-s", "--silence", help="静默模式，当把该值设为True时, 将不再在屏幕上输出日志信息，默认为False", default="False")
parser.add_argument("-m", "--mode", help="模式, 为sanger、tsanger、tsg中的一个， 默认为sanger", default="sanger")
# parser.add_argument("-p", "--port", help="端口号，通常为80或者2333，默认为80", default="80")
parser.add_argument("-f", "--fake", help="只发送文件路径，不传文件", default="False")
args = vars(parser.parse_args())

# if args["get_list"] =="T":
#     get_list(args["input_dir"],args["file_list"])
# else:
if not os.path.exists(args["input_dir"]):
    raise OSError("输入路径 {} 不存在".format(args["input_dir"]))

if os.path.isfile(args["input_dir"]):
    raise OSError("输入路径不能为文件，只能为文件夹")

if not os.path.exists(args["file_list"]):
    print("list文件{}不存在，使用输入目录{}生成list文件".format(args["file_list"], args["input_dir"]))
    get_list(args["input_dir"], args["file_list"])
    # raise OSError("列表文件 {} 不存在".format(args["file_list"]))

if args["identity_code"] == "none":
    raise ValueError("验证码不能为none")

# if not os.path.exists(args["file_list"]):
#     raise OSError("列表文件 {} 不存在".format(args["file_list"]))

if args["silence"] not in ["True", "False"]:
    raise ValueError("参数-s的值必须是True或者是False")

if args["mode"] not in ["sanger", "tsanger", "tsg", "sg"]:
    raise ValueError("参数-m的值必须是sanger或者是tsanger,tsg")

if args["silence"] == "False":
    stream_on = True
else:
    stream_on = False

# task = UploadTask(args["identity_code"],"", args["mode"], args["port"], stream_on)
task = UploadTask(os.path.abspath(args["input_dir"]), os.path.abspath(args["file_list"]), args["identity_code"], args["mode"], stream_on)
# info = task.get_task_info()
# task.get_file_info(os.path.abspath(args["input_dir"]),os.path.abspath(args["file_list"]))

if args["fake"] == "False":
    task.upload_files()

task.post_filesdata()
