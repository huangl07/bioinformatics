# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# modified by yuguo at 20170920

from __future__ import division
import urllib2
import urllib
# import httplib
import sys
import json
import requests
import os
# import pprint
import re
import random
import time
import hashlib
from requests_toolbelt import MultipartEncoder, MultipartEncoderMonitor
from basic import Basic


class FileLimiter(object):
    def __init__(self, file_obj, read_limit):
        self.read_limit = read_limit
        self.amount_seen = 0
        self.file_obj = file_obj
        self.len = read_limit

    def read(self, amount=-1):
        if self.amount_seen >= self.read_limit:
            return b''
        remaining_amount = self.read_limit - self.amount_seen
        data = self.file_obj.read(remaining_amount if amount < 0 else min(amount, remaining_amount))
        self.amount_seen += len(data)
        return data


class UploadTask(Basic):
    def __init__(self, input_dir, list_path, identity, mode, stream_on):
        super(UploadTask, self).__init__('', identity, mode, stream_on)
        self.logger.info(self.post_url)
        self.source_path = input_dir  # 完整路径 /mnt/.../sg-users/./updir
        tmp_list = re.split(r'[\\/]', self.source_path)
        # length = len(tmp_list)
        self.source_dir = tmp_list.pop()  # 输入文件夹 /updir

        # self._upload_url = self.get_upload_url(mode)  # 文件上传的php接口的地址
        self._file_info = list()
        self.target_path = self.post_verifydata()   # rerewrweset/files/m_190/10000001
        self.logger.info("--------source_path:" + self.source_path + "--source_dir" + self.source_dir)
        self.get_file_info(list_path)

        # self.no_prifix_path = ""  # rerewrweset/files/m_190/10000001
        # self.post_url = ""  # 与前端通信的的接口的地址，用于告诉前端文件已经上传， 以及具体上传的了何处

    # @property
    # def upload_url(self):
    #     return self._upload_url

    # def get_url(self, mode, port):
    #     self._url = "http://192.168.12.101:{}/app/dataexchange/upload_task".format(port)
    #     return self._url  # 前置的接口地址，用于验证和获取上传文件的路径

    # def get_upload_url(self, mode):
    #     self._upload_url = "http://192.168.12.101/upload.php"
    #     if mode == "tsg":
    #         self._upload_url = "http://192.168.12.101/tsgupload.php"
    #     return self._upload_url  # 接受文件的接口地址，php

    def get_file_info(self, list_path):
        # self.source_path = input_dir
        # tmp_list = re.split("/", self.source_path)

        with open(list_path, "rb") as r:
            # line = r.next().rstrip()
            # if not re.search("#source#", line):
            #     raise ValueError("输入的列表文件的格式错误！")
            # else:
            #     self.source_path = re.sub("#source#", "", line)  # updir
            line = r.next()
            for line in r:
                line = line.rstrip().split("\t")
                d = dict()
                d["path"] = line[0]
                d["size"] = re.sub("B", "", line[1])
                d["description"] = line[2]
                d["locked"] = line[3]
                rel_path = re.sub(self.source_path, "", d["path"]).lstrip("/")
                self.logger.info(rel_path)
                rel_path = re.sub(self.source_dir, "", rel_path).lstrip('/')
                # rel_path = os.path.join(soure_dir, rel_path).lstrip("/")
                d["rel_path"] = rel_path   # \t1.txt
                self.logger.info("--------source_path:" + self.source_path + "--source_dir:" + self.source_dir)
                self.logger.info("--------d[path]:" + d["path"] + "--d[rel_path]:" + d["rel_path"])
                d["target_path"] = os.path.join(self.config_path, self.target_path, self.source_dir, rel_path)  # /mnt/iuster/tsanger-data/rerewrweset/files/m_190/10000001/updir/t1.txt
                # self.logger.info("target: "+self.target_path+"confg_path:"+self.config_path+"source_dir"+self.source_dir)
                self.logger.info("---------d[target_path]:" + d["target_path"])

                self._file_info.append(d)
        return self._file_info
    # def get_post_url(self, mode):
    #     if self.post_url != "":
    #         return self.post_url
    #     if mode == "sanger":
    #         # self.post_url = "http://www.sanger.com/api/add_file_dir"
    #         self.post_url = "http://api.sanger.com/file/verify_filecode"
    #     if mode == "tsanger":
    #         # self.post_url = "http://www.tsanger.com/api/add_file_dir"
    #         self.post_url = "http://api.tsanger.com/file/verify_filecode"
    #     if mode == "tsg":
    #         self.post_url = "http://api.tsg.com/file/verify_filecode"
    #     return self.post_url

    # def get_task_info(self):
    #     """
    #     验证验证码， 获取上传的文件路径
    #     """
    #     data = urllib.urlencode({"ip": self.ip, "identity": self.identity, "user": self.user, "mode": self.mode})
    #     req = urllib2.Request(self.url, data)
    #     try:
    #         self.logger.info("用户: {}, 验证码: {}".format(self.user, self.identity))
    #         self.logger.info("正在与远程主机通信，获取项目信息")
    #         response = urllib2.urlopen(req)
    #     except (urllib2.HTTPError, urllib2.URLError, httplib.HTTPException) as e:
    #         self.logger.info(e)
    #         if self._port != "2333":
    #             try:
    #                 self.logger.info("尝试使用2333端口重新进行连接")
    #                 self.get_url(self.mode, "2333")
    #                 req = urllib2.Request(self.url, data)
    #                 response = urllib2.urlopen(req)
    #             except (urllib2.HTTPError, urllib2.URLError, httplib.HTTPException) as e:
    #                 self.logger.info(e)
    #                 sys.exit(1)
    #             else:
    #                 info = response.read()
    #         else:
    #             sys.exit(1)
    #     else:
    #         info = response.read()
    #
    #     info = json.loads(info)
    #     if not info["success"]:
    #         self.logger.info(info["info"])
    #         sys.exit(1)
    #     else:
    #         self.logger.info("通信成功，开始上传文件...")
    #         self._target_path = info["abs_path"]    # 经过config转化的绝对路径, /mnt/iluster/tsanger-data/+rel——path
    #         self.no_prifix_path = info["rel_path"]  # 即是数据库里的rel_dir的值， rerew...//
    #         return json.dumps(info)

    def post_verifydata(self):
        """
        请求url: http://api.sg.com/file/verify_filecode
        请求方式:post

        参数：
        $params = array(
            'code'     => 'QAEBXH|702692627ec061cf853b3317bfc1a776',
            'type'     => 'upload', //上传
            'dir_name' => 'abc',  //验证的文件夹名称
        );

        成功时,json结果:
        {
            "success": "true",
            "data": {
                "path": "rerewrweset/files/m_219/10008074"
            }
        }
        """
        my_data = dict()
        my_data['code'] = self.identity
        my_data['type'] = 'upload'
        my_data['dir_name'] = self.source_dir
        # my_data["data"] = dict()
        # my_data["data"]["files"] = list()
        # for d in self._file_info:
        #     my_data["data"]["files"].append({"path": d["submit_path"], "format": "", "description": d["description"], "locked": d["locked"], "size": d["size"]})
        # my_data["data"]["dirs"] = self._get_dir_struct()

        my_data["client"] = self.client
        my_data["nonce"] = str(random.randint(1000, 10000))
        my_data["timestamp"] = str(int(time.time()))
        x_list = [self.client_key, my_data["timestamp"], my_data["nonce"]]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        my_data["signature"] = sha1.hexdigest()
        # pprint.pprint(my_data)
        # my_data = json.dumps(my_data)
        request = urllib2.Request(self.post_url, urllib.urlencode(my_data))
        self.logger.info("与{}网站通信， 发送验证码和文件夹名称请求：{}".format(self.post_url, my_data))
        # self.logger.info("与sanger网站通信， 将上传结果传递至sanger网站上")
        self.logger.info("获取文件夹是否有重名，以及目标存放路径")
        try:
            response = urllib2.urlopen(request)
        except urllib2.HTTPError as e:
            self.logger.error(e)
            raise Exception(e)
        else:
            the_page = response.read()
            self.logger.info("Return Page:")
            self.logger.info(the_page)
            info = json.loads(the_page)
            self.logger.info(info)
            if info["success"] == 'false':
                self.logger.info(info["message"].encode('utf-8'))
                sys.exit(1)
            else:
                self.logger.info("通信成功，获得目标路径{}".format(info["data"]["path"]))
                target_path = info["data"]["path"]
                return target_path
                # self.logger.info(self.target_path)
            # my_return = json.loads(the_page)
            # if my_return["success"] in ["true", "True", True]:
            #     self.logger.info("文件上传已经全部结束！")
            # else:
            #     raise Exception("文件信息写入数据库失败：{}".format(my_return["message"]))

    def my_callback(self, monitor):
        upload_bite = monitor.bytes_read
        my_m = (upload_bite / 1024) / 1024
        print "已经上传{0:.2f}M".format(my_m)

    def upload_files(self):
        """
        上传文件
        """
        if len(self._file_info) == 0:
            raise ValueError("没有需要上传的文件")
        # multiple_files= [('sources', ("aa.txt", open(my_f, "rb"), "txt/plain")),
        #                 ("sources", ("bbb.txt", open(my_f2, "rb"), "txt/plain"))
        #    ]
        # 获取上传目标的文件夹名称， 为拼凑上传的目标路径做准备
        # tmp_list = re.split("/", self.source_path)
        # length = len(tmp_list)
        # soure_dir = tmp_list.pop()
        # for i in range(length):
        #     soure_dir = tmp_list.pop()
        #     if soure_dir == "":
        #         continue
        #     else:
        #         break
        # self.logger.info("upload: "+self.target_path)
        for d in self._file_info:
            # pprint.pprint(d)
            # rel_path = re.sub(self.source_path, "", d["path"]).lstrip("/")
            # self.logger.info(rel_path)
            # rel_path = re.sub(self.source_dir,"",rel_path).lstrip('/')
            # # rel_path = os.path.join(soure_dir, rel_path).lstrip("/")
            # d["rel_path"] = rel_path   # \t1.txt
            # self.logger.info("--------source_path:"+self.source_path+"--source_dir:"+self.source_dir)
            # self.logger.info("--------d[path]:"+d["path"]+"--d[rel_path]:"+d["rel_path"])
            # d["target_path"] = os.path.join(self.config_path, self.target_path, self.source_dir, rel_path)  # /mnt/iuster/tsanger-data/rerewrweset/files/m_190/10000001/updir/t1.txt
            # # self.logger.info("target: "+self.target_path+"confg_path:"+self.config_path+"source_dir"+self.source_dir)
            # self.logger.info("---------d[target_path]:"+d["target_path"])
            # psot_json = {"target_path": full_path, "mode": self.mode, "code": self.identity, "rel_dir": self.no_prifix_path}

            # my_file = {'sources': (os.path.basename(d["path"]), open(d["path"], "rb"), 'application/octet-stream'),
            #           'target': (None, json.dumps(psot_json), 'application/json')}
            # old version
            # m = MultipartEncoder(
            #     fields={'sources': (os.path.basename(d["path"]), open(d["path"], "rb"), 'application/octet-stream'),
            #             "target_path": full_path, "mode": self.mode, "code": self.identity, "rel_dir": self.no_prifix_path})

            # # fix version before change php file by yuguo
            # m = MultipartEncoder(
            #     fields={'sources': (os.path.basename(d["path"]), open(d["path"], "rb"), 'application/octet-stream'),
            #             "target_path": d["target_path"], "mode": self.mode, "code": self.identity, "rel_dir": self.target_path})

            # # fix version after change php file by yuguo

            m = MultipartEncoder(
                fields={'sources': (os.path.basename(d["path"]), open(d["path"], "rb"), 'application/octet-stream'),
                        "target_path": d["target_path"], "verify": "sanger-data-upanddown", "rel_dir": self.target_path + "/"})
            m = MultipartEncoderMonitor(m, self.my_callback)
            # d["target_path"] = full_path
            # if self.mode == "tsanger" or self.mode == "tsg":
            #     prefix = "tsanger:"
            # elif self.mode == "sanger":
            #     prefix = "sanger:"
            # d["submit_path"] = self.prefix_path + os.path.join(self.target_path, rel_path)  # 上传到硬盘的哪个位置
            self.logger.info("rel_dir:{}".format(self.target_path))
            self.logger.info("开始上传文件{}".format(d["path"]))
            # with open(d["path"], 'rb') as file_obj:
            #    upload = FileLimiter(file_obj, 40960)
            #    r = requests.post(self._upload_url, data=upload, headers={'Content-Type': 'application/octet-stream'}, json=psot_json)
            r = requests.post(self.upload_url, data=m, headers={'Content-Type': m.content_type})
            # print r.text
            if r.status_code == 200:
                self.logger.info("文件{}上传完成".format(d["path"]))
            else:
                print "ERROR status_code :" + r.text.encode('utf-8')
        # self.post_filesdata()

    def post_filesdata(self):
        """
        磁盘文件上传成功后，需要请求接口提供文件路径信息：

        请求url: http://api.sg.com/file/add_file_by_code
        请求方式:post
        参数：
        return array(
            'param' => array(
            'code' => 'QAEBXH|702692627ec061cf853b3317bfc1a776',
            'type' => 'upload',
            ),
            'base_path' => 'rerewrweset/files/m_219/10008074',  #前置路径
            'files'     => array(
            array(
                'path' => 'corr_network_calc/corr_network_centrality.txt', //路径
                'size' => '13750',   //大小
                'description' => 'OTU\\u5e8f\\u5217\\u7269\\u79cd\\u5206\\u7c7b\\u6587\\u4ef6', //描述
                'format' => 'taxon.seq_taxon',   //格式
            ),
            array(
                'path' => 'corr_network_calc/corr_network_by_cut.txt',
                'size' => '303363',
                'description' => 'OTU\u5e8f\u5217\u7269\u79cd\u5206\u7c7b\u6587\u4ef6',
                'format' => 'txt',
            ),
            array(
                'path' => 'abc.txt',
                'size' => '1000',
                'description' => 'rererwew',
                'format'      => 'txt',
            )
            ),
            //给相关目录添加描述信息
            'dirs' => array(
            array(
                'path' => 'rerewrweset/files/m_219/10008074/tsanger_3/report_results', //路径
                'description' => '\u57fa\u7840\u5206\u6790\u7ed3\u679c\u6587\u4ef6\u5939', //描述
            ),
            array(
                'path' => 'corr_network_calc',
                'description' => '物种相关性网络分析结果输出目录',
            ),
            ),
        );
        """
        self.logger.info("post_data: " + self.target_path)
        post_data = dict()
        post_data['param'] = {'code': self.identity, 'type': 'upload'}
        post_data['base_path'] = self.prefix_path + os.path.join(self.target_path, self.source_dir)
        post_data["files"] = list()
        post_data["dirs"] = list()
        # my_data["data"]["files"] = list()
        # for d in self._file_info:
        #     my_data["data"]["files"].append({"path": d["submit_path"], "format": "", "description": d["description"], "locked": d["locked"], "size": d["size"]})
        # my_data["data"]["dirs"] = self._get_dir_struct()
        for d in self._file_info:
            post_data["files"].append({"path": d["rel_path"], "format": "", "description": d["description"], "lock": d["locked"], "size": d["size"]})
        # my_data["dirs"] = self._get_dir_struct()

        my_data = dict()
        my_data["client"] = self.client
        my_data["nonce"] = str(random.randint(1000, 10000))
        my_data["timestamp"] = str(int(time.time()))
        x_list = [self.client_key, my_data["timestamp"], my_data["nonce"]]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        my_data["signature"] = sha1.hexdigest()
        # pprint.pprint(my_data)
        my_data["sync_files"] = json.dumps(post_data)
        request = urllib2.Request(self.post_add_url, urllib.urlencode(my_data))
        self.logger.info("与{}网站通信， 发送验证码和文件列表：{}".format(self.post_add_url, my_data))
        try:
            response = urllib2.urlopen(request)
        except urllib2.HTTPError as e:
            self.logger.error(e)
            raise Exception(e)
        else:
            the_page = response.read()
            self.logger.info("Return Page:")
            self.logger.info(the_page)
            my_return = json.loads(the_page)
            if my_return["success"] in ["true", "True", True]:
                self.logger.info("文件上传已经全部结束！")
            else:
                raise Exception("发送网站文件信息失败：{}".format(my_return["message"].encode('utf-8')))

    # def _get_dir_struct(self):
    #     my_dict = list()
    #     dir_list = list()
    #     for d in self._file_info:
    #         dir_list.append(os.path.dirname(d["rel_path"]))
    #     dir_list = list(set(dir_list))
    #     if self.mode == "tsanger"  or self.mode == "tsg":
    #         prefix = "tsanger:"
    #     elif self.mode == "sanger":
    #         prefix = "sanger:"
    #     for l in dir_list:
    #         my_dict.append({"path": prefix + os.path.join(self.no_prifix_path, l), "format": "", "description": "", "locked": "", "size": ""})
    #     return my_dict
