# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# modified by yuguo at 20170719
# modified by yuguo at 20170920

from __future__ import division
import urllib2
import urllib
import httplib
import sys
import json
import requests
import os
import errno
import random
import time
import hashlib
from basic import Basic


class DownloadTask(Basic):
    """
    用于内部下载平台上的文件
    """
    def __init__(self, identity, outpath, mode, stream_on):
        super(DownloadTask, self).__init__(outpath,identity, mode, stream_on)
        # self.post_url = ""  # 与前端通信的的接口的地址
        # self._download_url = self.get_download_url(mode) # 存放数据服务器地址
        self._file_list = list()
        self.outpath = outpath

    # @property
    # def download_url(self):
    #     return self._download_url

    # def get_url(self, mode, port):
    #     if mode == "sanger":
    #         self._url = "http://192.168.12.101:{}/app/dataexchange/download_task".format(port)
    #     elif mode == "tsanger" or mode == "tsg":
    #         self._url = "http://192.168.12.102:{}/app/dataexchange/download_task".format(port)
    #     return self._url

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

    def get_download_url(self, mode):
        if mode == "sanger":
            self._download_url = "http://192.168.12.101/download.php"
        if mode == "tsanger" or mode == "tsg":
            self._download_url = "http://192.168.12.101/tsgdownload.php"
        return self._download_url

    def post_verifydata(self):
        """
        验证下载码：
        请求url: http://api.sg.com/file/verify_filecode
        请求方式:post
        参数：
        $params = array(
            'code'     => 'ATWKZN|e0cac412c4956c0879f2025b51d2024b',
            'type'     => 'download', //下载
        );
        成功时返回的数据：
        {
            "success": "true",
            "data": {
                "files": [
                    {
                        "file_name": "seqs_tax_assignments.txt",
                        "current_path": "sanger_15004/cmd_112_1493109842/output/Tax_assign/",
                        "disk_url": "rerewrweset/files/m_193/10007924/i-sanger_14621/cmd_112_1493109842/output/Tax_assign/seqs_tax_assignments.txt"
                    },
                    {
                        "file_name": "valid_sequence.txt",
                        "current_path": "sanger_15004/cmd_112_1493109842/output/QC_stat/",
                        "disk_url": "rerewrweset/files/m_193/10007924/i-sanger_14621/cmd_112_1493109842/output/QC_stat/valid_sequence.txt"
                    },
                    ...
                ]
            }
        }
        """
        my_data = dict()
        my_data['code'] = self.identity
        my_data['type'] = 'download'
        # my_data["data"] = dict()
        # my_data["data"]["files"] = list()
        # for d in self._file_info:
        #     my_data["data"]["files"].append({"path": d["submit_path"], "format": "", "description": d["description"], "locked": d["locked"], "size": d["size"]})
        # my_data["data"]["dirs"] = self._get_dir_struct()
        my_data['client'] = self.client
        my_data['nonce'] = str(random.randint(1000, 10000))
        my_data['timestamp'] = str(int(time.time()))
        x_list = [self.client_key, my_data['timestamp'], my_data['nonce']]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        my_data['signature'] = sha1.hexdigest()
        request = urllib2.Request(self.post_url, urllib.urlencode(my_data))
        self.logger.info("与{}网站通信， 发送下载验证请求：{}".format(self.post_url, my_data))
        try:
            response = urllib2.urlopen(request)
        except urllib2.HTTPError as e:
            self.logger.error(e)
            raise Exception(e)
        else:
            the_page = response.read()
            self.logger.info("Return Page:")
            self.logger.info(the_page)
            # my_return = json.loads(the_page)
            # print my_return
            info = json.loads(the_page)
            if not info["success"]:
                self.logger.info(info["message"])
                sys.exit(1)
            else:
                self.logger.info("通信成功，获得文件列表")
                self._file_list = info["data"]["files"]
                print info["data"]

            # if my_return["success"] in ["true", "True", True]:
            #     self.logger.info("文件上传已经全部结束！")
            # else:
            #     raise Exception("文件信息写入数据库失败：{}".format(my_return["message"]))



    def download_files(self):
        """
        遍历文件列表， 并下载文件
            "files": [
        {
            "file_name": "seqs_tax_assignments.txt",
            "current_path": "sanger_15004/cmd_112_1493109842/output/Tax_assign/",
            "disk_url": "rerewrweset/files/m_193/10007924/i-sanger_14621/cmd_112_1493109842/output/Tax_assign/seqs_tax_assignments.txt"
        }
        """
        total_sum = len(self._file_list)
        count = 1
        for f_info in self._file_list:
            file_name = os.path.basename(f_info["file_name"])
            dir_name = os.path.dirname(f_info["current_path"])
            local_dir = os.path.join(self.outpath, dir_name)
            local_file = os.path.join(local_dir, file_name)
            disk_url = self.config_path+f_info['disk_url']
            # biocluster接口：file_list.append([full_path, file_size, rel_path, 2])
            # file_name = os.path.basename(f_info[0])
            # dir_name = os.path.dirname(f_info[2])
            # local_dir = os.path.join(self.outpath, dir_name)
            # local_file = os.path.join(local_dir, file_name)

            # self.logger.info("正在下载第 {}/{} 个文件: {}, 文件大小{}".format(count, total_sum, file_name, f_info[1]))
            self.logger.info("正在下载第 {}/{} 个文件: {}...".format(count, total_sum, file_name))
            count += 1
            #post_info = urllib.urlencode({'indentity_code': self.identity, 'file': f_info[0], 'mode': self.mode})
            post_info = urllib.urlencode({'file': disk_url,'verify': "sanger-data-upanddown"})
            request = urllib2.Request(self.download_url, post_info)
            self.logger.info("request:{}".format(self.download_url))
            try:
                u = urllib2.urlopen(request)
                self.logger.info("{}".format(u))
                try:
                    os.makedirs(local_dir)
                except OSError as exc:
                    if exc.errno != errno.EEXIST:
                        raise exc
                    else:
                        pass
            except (urllib2.HTTPError, urllib2.URLError, httplib.HTTPException) as e:
                self.logger.info(e)
                continue
            meta = u.info()
            file_size = int(meta.getheaders("Content-Length")[0])
            file_size_dl = 0
            block_sz = 51200
            f = open(local_file, "wb")
            while True:
                buffer = u.read(block_sz)
                if not buffer:
                    break

                file_size_dl += len(buffer)
                f.write(buffer)
                status = "{:10d}  [{:.2f}%]".format(file_size_dl, file_size_dl * 100 / file_size)
                self.logger.info(status)
            f.close()

    # def get_task_info(self):
    #     """
    #     获取编码相对应的task_id, 以及相关的目录结构
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

    #     info = json.loads(info)
    #     if not info["success"]:
    #         self.logger.info(info["info"])
    #         sys.exit(1)
    #     else:
    #         self.logger.info("通信成功，任务id为：{}。开始下载文件...".format(info["task_id"]))
    #         self._file_list = info["data"]
    #         print info["data"]
