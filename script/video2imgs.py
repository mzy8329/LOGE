import cv2

def video_to_images(video_path, image_path):
  cap = cv2.VideoCapture(video_path)
  fps = cap.get(cv2.CAP_PROP_FPS)
  frame_count = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
  count = 0
  while(cap.isOpened()):
      ret, frame = cap.read()
      if ret == True:
          cv2.imwrite(image_path + "/frame_" + str(count) + ".jpg", frame)
          count += 1
      else:
          break

  cap.release()
  cv2.destroyAllWindows()

video_path = "data/videos/underwater.mp4"
image_path = "data/imgs"
video_to_images(video_path, image_path)